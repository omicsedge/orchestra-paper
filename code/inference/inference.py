import json
import logging
import os
import tempfile
from pathlib import Path

import click
import dask.dataframe as dd
import pandas as pd

from .app import Transformer, convert_to_ancestry_format, predict_fn, vcf_preprocess


def run_command(cmd: str) -> None:
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise RuntimeError(f"Command failed with exit code {exit_code}: {cmd}")


def weight_to_str(x):
    x = x / x.sum()
    x = x.sort_values(ascending=False)
    return ", ".join(f"{idx[2]}: {v*100:0.2f}" for idx, v in x.items())


@click.command()
@click.option("-panel", "-p", type=Path, help="Inference panel path")
@click.option(
    "--output-dir",
    "-o",
    type=Path,
    help="Output directory",
    default="/results/inference",
)
@click.option("--model-path", "-m", type=Path, help="Model path")
def inference(panel: Path, output_dir: Path, model_path: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        run_command(
            f"bash /code/inference/scripts/preprocessing.sh {panel} {temp_path}"
        )

        submodels = model_path / "sub-models"

        for submodel in submodels.iterdir():
            logging.info(f"Running submodel {submodel}")
            logging.info(f"Output dir: {output_dir}")

            model_dir = submodel / "base" / "model.tar.gz"
            smodel_dir = submodel / "smooth" / "model.tar.gz"

            (temp_path / submodel.name).mkdir(parents=True, exist_ok=True)
            run_command(f"tar -xzf {model_dir} -C {temp_path / submodel.name}")
            run_command(f"tar -xzf {smodel_dir} -C {temp_path / submodel.name}")
            model_dir = temp_path / submodel.name / "model"

            sample_list_dir = output_dir / "window_results" / submodel.name
            vcf_preprocess(
                input_dir=temp_path / "dataset",
                model_dir=model_dir,
                output_dir=sample_list_dir,
            )
            processor = Transformer(model_dir)

            for sample_list_path in sample_list_dir.iterdir():
                sample_list_file = sample_list_path.name

                data = pd.read_csv(sample_list_path, sep="\t", header=None)

                # predict
                logging.info(f"Running base layer on {sample_list_file}")
                predictions = processor.transform(data)
                (output_dir / "base_results" / submodel.name).mkdir(
                    parents=True, exist_ok=True
                )
                predictions.to_csv(
                    output_dir / "base_results" / submodel.name / sample_list_file,
                    header=False,
                    index=False,
                    sep="\t",
                )

                # predict smoothing layer
                y_pred = predict_fn(predictions, model_dir)
                (output_dir / "smooth_results" / submodel.name).mkdir(
                    parents=True, exist_ok=True
                )
                y_pred.to_csv(
                    output_dir / "smooth_results" / submodel.name / sample_list_file,
                    header=False,
                    index=False,
                    sep="\t",
                )

            df_base = convert_to_ancestry_format(
                output_dir / "base_results" / submodel.name,
                model_path / "artifacts" / "population_map.tsv",
                submodel / "artifacts" / "parameters.json",
                pred_by_argmin=True,
            )
            (output_dir / "summary_results").mkdir(parents=True, exist_ok=True)
            df_base.to_csv(
                output_dir / "summary_results" / f"base_samples.{submodel.name}.tsv.gz",
                sep="\t",
                index=False,
            )
            df_smooth = convert_to_ancestry_format(
                output_dir / "smooth_results" / submodel.name,
                model_path / "artifacts" / "population_map.tsv",
                submodel / "artifacts" / "parameters.json",
                pred_by_argmin=False,
            )
            df_smooth.to_csv(
                output_dir
                / "summary_results"
                / f"smooth_samples.{submodel.name}.tsv.gz",
                sep="\t",
                index=False,
            )

        df_base = dd.read_csv(
            output_dir / "summary_results" / "base_samples.*.tsv.gz",
            sep="\t",
            blocksize=None,
            dtype={"sample_id": "object"},
        ).compute()
        df_base["layer"] = "base"

        df_smooth = dd.read_csv(
            output_dir / "summary_results" / "smooth_samples.*.tsv.gz",
            sep="\t",
            blocksize=None,
            dtype={"sample_id": "object"},
        ).compute()
        df_smooth["layer"] = "smooth"

        # Performing basic analysis of the results.
        df_base["layer"] = "base"
        df_smooth["layer"] = "smooth"

        df = pd.concat([df_smooth, df_base])
        del df_base, df_smooth

        with open("/code/inference/chr_map.json", "rt") as fin:
            chr_fractions = json.load(fin)

        chr_fractions = {int(k): v for k, v in chr_fractions.items()}
        window_lengths = df.groupby(["chrom"])["window"].nunique().to_dict()
        chr_map = {
            k: v / window_lengths[k]
            for k, v in chr_fractions.items()
            if k in window_lengths
        }

        df["weight"] = df.chrom.map(chr_map)

        df_ = (
            df.groupby(["sample_id", "layer", "pred"])
            .agg({"weight": "sum"})
            .groupby(["sample_id", "layer"])
            .agg({"weight": weight_to_str})
            .rename(columns={"weight": "predicted ancestry"})
            .reset_index()
            .sort_values(["sample_id", "layer"])
        )
        df_.to_csv(
            output_dir / "summary_results" / f"ancestry.tsv",
            sep="\t",
            index=False,
        )


if __name__ == "__main__":
    inference()
