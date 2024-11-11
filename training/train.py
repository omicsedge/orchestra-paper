"""
Script for model training
"""

import logging
import os
from pathlib import Path
from typing import List, Tuple

import click
from joblib import Parallel, delayed
from scripts import (
    AncestryModel,
    evaluate,
    merge_base_layer,
    predict_window,
    smooth_inference,
    train_base_layer,
    train_smooth_layer,
    zarr_to_parquet,
)


@click.command()
@click.option(
    "--simulated-data-path",
    "-sd",
    help="Simulated dataset path",
    type=Path,
)
@click.option("--start-chr", "-sc", type=int, help="Start chromosome")
@click.option("--end-chr", "-ec", type=int, help="End chromosome")
@click.option(
    "--output-dir",
    "-o",
    help="Output for data",
    type=Path,
)
@click.option(
    "--window-size",
    "-ws",
    help="Window size",
    type=int,
    default=600,
)
@click.option(
    "--level",
    "-l",
    help="Level",
    type=int,
    default=3,
)
@click.option(
    "--version",
    "-v",
    help="Version",
    type=str,
    default="0.01",
)
@click.option(
    "--cohorts",
    "-c",
    multiple=True,
    default=["train", "test"],
    show_default=True,
    help="List of cohorts to process.",
)
def train(
    simulated_data_path: Path,
    start_chr: int,
    end_chr: int,
    output_dir: Path,
    window_size: int,
    level: int,
    version: str,
    cohorts: List[str],
):
    logging.info("Run preprocessing, converting zarr to parquet + windows")
    logging.info(f"Chromosomes: {start_chr}-{end_chr}")

    output_dir = output_dir / f"m{start_chr}.{end_chr}"

    output_dir.mkdir(parents=True, exist_ok=True)
    chromosomes: List[int] = list(range(start_chr, end_chr + 1))
    cpu_count: int = os.cpu_count()
    cohorts: Tuple[str, str] = ("train", "test")
    dataset: Path = output_dir / f"w{window_size}_parquet"

    # ConvertZarrToParquet
    Parallel(n_jobs=cpu_count, prefer="threads")(
        delayed(zarr_to_parquet)(
            pop_map=str(simulated_data_path / "population_map.tsv"),
            input_dir=simulated_data_path / cohort / "zarr-files" / f"chr{chrom}",
            output_dir=dataset / cohort,
            # Ensures separate directories for each cohort
            chromosome=chrom,
            window_size=window_size,
            level=level,
        )
        for cohort in cohorts
        for chrom in chromosomes
        # This ensures that every chromosome for each cohort is processed
    )

    # TrainBaseModel: ust create a output.tar.gz from windows info and model.tar.gz from gen0 of previous step.
    train_base_layer(
        input_dir=dataset / "train", output_dir=output_dir / "base_layer_model"
    )

    # PredictBaseLayer:
    Parallel(n_jobs=max(1, cpu_count - 2), prefer="threads")(
        delayed(predict_window)(
            data_file=data_file,
            base_layer_model=output_dir / "base_layer_model" / "model" / f"chr{chrom}",
            output_dir=output_dir / "base_inference" / cohort / f"chr{chrom}",
            training_set=cohort == "train",
        )
        for cohort in cohorts
        for chrom in chromosomes
        for data_file in (dataset / cohort / f"chr{chrom}").rglob("data.parquet")
    )

    # MergeBasePrediction
    for cohort in cohorts:
        merge_base_layer(
            input_dir=output_dir / "base_inference" / cohort,
            output_dir=output_dir / "base_inference" / cohort,
        )

    # TrainSmoothingModel
    train_smooth_layer(
        pop_map=simulated_data_path / "population_map.tsv",
        training_dir=output_dir / "base_inference" / "train",
        validation_dir=output_dir / "base_inference" / "test",
        output_dir=output_dir / "smooth_layer_model",
        epochs=100,
        level=level,
    )
    # PredictSmoothLayer-test
    smooth_inference(
        smooth_layer_model=output_dir / "smooth_layer_model",
        base_layer_pred=output_dir / "base_inference" / "test" / "predictions.tsv.gz",
        base_layer_params=output_dir / "base_inference" / "test" / "parameters.json",
        output_dir=output_dir / "smooth_inference_test",
    )
    # SaveModel
    AncestryModel(
        population_map_file=simulated_data_path / "population_map.tsv",
        model_path=output_dir.parent / version,
        sub_model_name="chr" + ".".join(map(str, chromosomes)),
        windows_info_file=output_dir / "base_layer_model" / "windows_info.tsv",
        simulated_params_file=output_dir / "smooth_layer_model" / "parameters.json",
        base_model_path=output_dir / "base_layer_model" / "model",
        smooth_model_path=output_dir / "smooth_layer_model" / "model.pt",
    ).save()

    # Evaluate-base-test
    evaluate(
        input_dir=output_dir / "base_inference" / "test",
        output_dir=output_dir / "evaluation" / "base",
        is_recomb=True,
    )
    # Evaluate-final-test
    evaluate(
        input_dir=output_dir / "smooth_inference_test",
        output_dir=output_dir / "evaluation" / "smooth",
    )


if __name__ == "__main__":
    train()
