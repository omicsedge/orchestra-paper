import logging
import os
import tempfile
from pathlib import Path

import click
import numpy as np

# import dask.dataframe as dd
import pandas as pd
from scripts import Transformer, vcf_preprocess


def run_command(cmd: str) -> None:
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise RuntimeError(f"Command failed with exit code {exit_code}: {cmd}")


@click.command()
@click.option("-panel", "-p", type=Path, help="Inference panel path")
@click.option("--name", "-n", type=str, help="Name")
@click.option("--output-dir", "-o", type=Path, help="Output directory")
def inference(panel: Path, name: str, output_dir: Path):
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        run_command(f"bash scripts/preprocessing.sh {panel} {temp_path}")

        training_path = Path("/output_training")

        for submodel in training_path.iterdir():
            logging.info(f"Running submodel {submodel}")
            logging.info(f"Output dir: {output_dir}")

            sample_list_dir = output_dir / "window_results" / submodel.parts[-1]
            model_dir = submodel / "base_layer_model" / "model"
            vcf_preprocess(
                input_dir=temp_path / "dataset",
                model_dir=model_dir,
                output_dir=sample_list_dir,
            )
            processor = Transformer(model_dir)

            for sample_list_file in sample_list_dir.iterdir():
                logging.info(f"Processing {sample_list_file}")

                data = pd.read_csv(sample_list_file, sep="\t", header=None)

                # predict
                processor.transform(data, output_file=f"{sample_list_file}.out")

                y_pred = predict_fn(X, model_dir)



    # jobs = []
    # for SUBMODEL_NAME in SUBMODELS:
    #     print(f"Running submodel {SUBMODEL_NAME}")

    #     job_name = f"smooth-{uuid.uuid4()}"
    #     jobs.append(job_name)

    #     # get directory
    #     SMOOTH_OUTPUT_DIR = f"{OUTPUT_DIR.smooth_results}{SUBMODEL_NAME}/"
    #     print(SMOOTH_OUTPUT_DIR)

    #     # create model
    #     model_registry = get_model_registry(ANCESTRY_MODEL, [SUBMODEL_NAME])

    #     model = PyTorchModel(
    #         model_data=model_registry["sub-models"][SUBMODEL_NAME]["model"][
    #             "smooth_model_uri"
    #         ],
    #         role=ROLE,
    #         framework_version="1.10",
    #         py_version="py38",
    #         source_dir="src/",
    #         entry_point="inference_user.py",
    #     )

    #     # batch transformer
    #     transformer = model.transformer(
    #         output_path=SMOOTH_OUTPUT_DIR,
    #         instance_count=1,
    #         instance_type="ml.m5.large",
    #         accept="text/tsv",
    #         strategy="MultiRecord",
    #         max_payload=4,
    #     )

    #     # run transformer
    #     BASE_OUTPUT_DIR = f"{OUTPUT_DIR.base_results}{SUBMODEL_NAME}/"
    #     transformer.transform(
    #         data=BASE_OUTPUT_DIR,
    #         content_type="text/tsv",
    #         split_type="Line",
    #         job_name=job_name,
    #         model_client_config={"InvocationsMaxRetries": 0},
    #         logs=False,
    #         wait=False,
    #     )

    # client = boto3.client("sagemaker")
    # errors = []

    # while len(jobs):
    #     print(f"{len(jobs)} left")
    #     for job_name in jobs:
    #         response = client.describe_transform_job(
    #             TransformJobName=job_name,
    #         )["TransformJobStatus"]

    #         if response in ["InProgress", "Stopping"]:
    #             time.sleep(60)
    #             continue

    #         if response in ["Failed", "Stopped"]:
    #             errors.append(job_name)

    #         indx = jobs.index(job_name)
    #         jobs.pop(indx)
    #         break

    # for error in errors:
    #     raise Exception(f"Smooth transform job {error} failed.")

    # print("Next step")

    # # # Base & smoothing layers post-processing

    # # convert base layer results to summary format
    # # generate summary of the smoothing layer output
    # for SUBMODEL_NAME in SUBMODELS:
    #     print(f"Running submodel {SUBMODEL_NAME}")

    #     # model registry
    #     model_registry = get_model_registry(ANCESTRY_MODEL, [SUBMODEL_NAME])

    #     parameters_file = model_registry["sub-models"][SUBMODEL_NAME]["artifacts"][
    #         "parameters"
    #     ]

    #     os.system("rm -rf temp/")
    #     os.system("mkdir -p temp/")
    #     os.system(f"aws s3 cp {parameters_file} temp/ --quiet")

    #     with open(f"temp/parameters.json", "rt") as fin:
    #         params = json.load(fin)

    #     for layer in ["base", "smooth"]:
    #         # base output path
    #         path = getattr(OUTPUT_DIR, f"{layer}_results")
    #         LAYER_OUTPUT_DIR = f"{path}{SUBMODEL_NAME}/"

    #         # read data and convert to regular format
    #         os.system(
    #             f"aws s3 cp --recursive {LAYER_OUTPUT_DIR} temp/{layer}/intermediate/ --quiet"
    #         )

    #         # read base results & save
    #         df_base = convert_to_ancestry_format(
    #             dd.read_csv(
    #                 f"temp/{layer}/intermediate/*",
    #                 sep="\t",
    #                 header=None,
    #                 dtype={2: "object"},
    #             ).compute(),
    #             model_registry["artifacts"]["population_map_uri"],
    #             params,
    #             pred_by_argmin=bool(layer == "base"),
    #         )
    #         df_base.to_csv(
    #             f"{OUTPUT_DIR.base}{OUTPUT_DIR.path}summary-results/{layer}_samples.{SUBMODEL_NAME}.tsv.gz",
    #             sep="\t",
    #             index=False,
    #         )
    #         del df_base

    # # ### Performing basic analysis of the results (saved to ancestry.tsv)

    # os.system("mkdir -p temp/")
    # os.system(
    #     f"aws s3 cp {OUTPUT_DIR.base}{OUTPUT_DIR.path}summary-results/ temp/ --recursive"
    # )

    # df_base = dd.read_csv(
    #     "temp/base_samples.*.tsv.gz",
    #     sep="\t",
    #     blocksize=None,
    #     dtype={"sample_id": "object"},
    # ).compute()
    # df_base["layer"] = "base"

    # df_smooth = dd.read_csv(
    #     "temp/smooth_samples.*.tsv.gz",
    #     sep="\t",
    #     blocksize=None,
    #     dtype={"sample_id": "object"},
    # ).compute()
    # df_smooth["layer"] = "smooth"

    # df = pd.concat([df_smooth, df_base])
    # del df_base, df_smooth

    # with open("chr_map.json", "rt") as fin:
    #     chr_fractions = json.load(fin)

    # chr_fractions = {int(k): v for k, v in chr_fractions.items()}
    # window_lengths = df.groupby(["chrom"])["window"].nunique().to_dict()
    # chr_map = {
    #     k: v / window_lengths[k]
    #     for k, v in chr_fractions.items()
    #     if k in window_lengths
    # }

    # df["weight"] = df.chrom.map(chr_map)

    # def weight_to_str(x):
    #     x = x / x.sum()
    #     x = x.sort_values(ascending=False)
    #     return ", ".join(f"{idx[2]}: {v*100:0.2f}" for idx, v in x.items())

    # df_ = (
    #     df.groupby(["sample_id", "layer", "pred"])
    #     .agg({"weight": "sum"})
    #     .groupby(["sample_id", "layer"])
    #     .agg({"weight": weight_to_str})
    #     .rename(columns={"weight": "predicted ancestry"})
    #     .reset_index()
    #     .sort_values(["sample_id", "layer"])
    # )
    # df_.to_csv(f"{NAME}-ancestry.tsv", sep="\t", index=False)

    # os.system("rm -rf temp/")


if __name__ == "__main__":
    inference()
