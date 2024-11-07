import json
import math
import os
import time
import uuid
from urllib.parse import urlparse

import boto3
import dask.dataframe as dd
import pandas as pd
from sagemaker.model import Model
from sagemaker.pytorch import PyTorchModel
from scripts import (
    AncestryModel,
    convert_to_ancestry_format,
    get_model_registry,
    preprocess,
)

# # Configuration

# ## Environment


BUCKET = "rnd-sandbox-datasets"
BASE_IMAGE_URI = "778090103881.dkr.ecr.us-east-1.amazonaws.com/sagemaker-training-containers/ancestry-recomb:latest"
ROLE = "arn:aws:iam::778090103881:role/ancestry_sagemaker_execution"

# Input CSV with location of VCF files.
NAME = "east-indonesians-alt7"
SAMPLES_SET_PATH = f"s3://rnd-sandbox-datasets/inference/cohorts/east-indonesians/"
OUTPUT_DIR_PATH = f"s3://rnd-sandbox-datasets/inference/results/{NAME}/"
MODEL_DIR_PATH = (
    f"s3://rnd-sandbox-datasets/sagemaker/model-registry/lai/combined/{NAME}/"
)
MODEL_VERSION = "0.01"


class S3Path:
    def __init__(self, url: str) -> None:
        parsed_url = urlparse(url)
        self.url = url
        self.path = parsed_url.path[1:]
        self.bucket = parsed_url.hostname
        self.sagemaker_registry_path = "sagemaker/model-registry/lai/"

    @property
    def input_files(self) -> str:
        return os.path.join(self.path, "input-files/")

    @property
    def window_results(self) -> str:
        return os.path.join(self.base, self.path, "window-results/")

    @property
    def base_results(self) -> str:
        return os.path.join(self.base, self.path, "base-results/")

    @property
    def smooth_results(self) -> str:
        return os.path.join(self.base, self.path, "smooth-results/")

    @property
    def model_name(self) -> str:
        return self.path.lstrip(self.sagemaker_registry_path).strip("/")

    @property
    def model_base_path(self) -> str:
        return f"{self.base}{self.sagemaker_registry_path}"

    @property
    def base(self) -> str:
        return f"s3://{self.bucket}/"


SAMPLES_SET = S3Path(SAMPLES_SET_PATH)
OUTPUT_DIR = S3Path(OUTPUT_DIR_PATH)
MODEL_DIR = S3Path(MODEL_DIR_PATH)


def write_file(data, idx):
    path = os.path.join(OUTPUT_DIR.input_files, f"input-data-{idx}.csv")
    client.put_object(Body=data, Bucket=OUTPUT_DIR.bucket, Key=path)
    print(f"upload: s3://{OUTPUT_DIR.bucket}/{path}")


client = boto3.client("s3")
result = client.list_objects(
    Bucket=SAMPLES_SET.bucket,
    Prefix=SAMPLES_SET.path,
    Delimiter="/",
)

folders = result["CommonPrefixes"]
lines_per_file = math.ceil(len(folders) / 5)
counter = 1
lines = ""
for i, o in enumerate(folders):
    lines += (
        ",".join([o["Prefix"].split("/")[-2], SAMPLES_SET.bucket, o["Prefix"]]) + "\n"
    )
    if (i + 1) % lines_per_file == 0:
        write_file(lines, counter)
        counter += 1
        lines = ""

if lines:
    write_file(lines, counter)

# ## Read the Model

ANCESTRY_MODEL = AncestryModel(
    MODEL_DIR.model_base_path,
    MODEL_DIR.model_name,
    MODEL_VERSION,
)
SUBMODELS = ANCESTRY_MODEL._submodels

# # Preprocess
#
# Convert dataset into windows.


jobs = []
for SUBMODEL_NAME in SUBMODELS:
    print(f"Running submodel {SUBMODEL_NAME}")

    # base output path
    WINDOW_OUTPUT_DIR = f"{OUTPUT_DIR.window_results}{SUBMODEL_NAME}/"
    print(WINDOW_OUTPUT_DIR)

    # model
    model_registry = get_model_registry(ANCESTRY_MODEL, [SUBMODEL_NAME])
    base_job_name = "vcfwindowizer"
    job_name = f"{base_job_name}-{uuid.uuid4()}"
    jobs.append(job_name)

    preprocess(
        input_dir=f"{OUTPUT_DIR.base}{OUTPUT_DIR.input_files}",
        model_dir=model_registry["sub-models"][SUBMODEL_NAME]["model"][
            "base_model_uri"
        ],
        output_dir=WINDOW_OUTPUT_DIR,
    )

# ### Waiting for results from preprocessing step


client = boto3.client("sagemaker")
errors = []

while len(jobs):
    print(f"{len(jobs)} left")
    for job_name in jobs:
        response = client.describe_processing_job(
            ProcessingJobName=job_name,
        )["ProcessingJobStatus"]

        if response in ["InProgress", "Stopping"]:
            time.sleep(60)
            continue

        if response in ["Failed", "Stopped"]:
            errors.append(job_name)

        indx = jobs.index(job_name)
        jobs.pop(indx)
        break

for error in errors:
    raise Exception(f"Processing job {error} failed.")

print("Next step")


# # Base Layer Inference


jobs = []
for SUBMODEL_NAME in SUBMODELS:
    print(f"Running submodel {SUBMODEL_NAME}")
    job_name = f"base-{uuid.uuid4()}"
    jobs.append(job_name)

    # base output path
    BASE_OUTPUT_DIR = f"{OUTPUT_DIR.base_results}{SUBMODEL_NAME}/"
    print(BASE_OUTPUT_DIR)

    # create model
    model_registry = get_model_registry(ANCESTRY_MODEL, [SUBMODEL_NAME])
    model = Model(
        image_uri=BASE_IMAGE_URI,
        model_data=model_registry["sub-models"][SUBMODEL_NAME]["model"][
            "base_model_uri"
        ],
        role=ROLE,
    )

    # create transformer
    transformer = model.transformer(
        instance_type="ml.m5.12xlarge",  # NOTE: for more than 2 chromosomes use ml.m5.2xlarge
        instance_count=5,
        output_path=BASE_OUTPUT_DIR,
        strategy="MultiRecord",
        max_payload=24,  # NOTE: for ml.m5.2xlarge use 4 as the payload
        max_concurrent_transforms=1,
        env={"MODEL_SERVER_TIMEOUT": "3600"},
    )

    # transform data
    WINDOW_OUTPUT_DIR = f"{OUTPUT_DIR.window_results}{SUBMODEL_NAME}/"
    transformer.transform(
        data=WINDOW_OUTPUT_DIR,
        content_type="text/csv",
        split_type="Line",
        model_client_config={
            "InvocationsMaxRetries": 0,
            "InvocationsTimeoutInSeconds": 3600,
        },
        job_name=job_name,
        logs=False,
        wait=False,
    )


client = boto3.client("sagemaker")
errors = []
while len(jobs):
    print(f"{len(jobs)} left")
    for job_name in jobs:
        response = client.describe_transform_job(
            TransformJobName=job_name,
        )["TransformJobStatus"]

        if response in ["InProgress", "Stopping"]:
            time.sleep(60)
            continue

        if response in ["Failed", "Stopped"]:
            errors.append(job_name)

        indx = jobs.index(job_name)
        jobs.pop(indx)
        break

for error in errors:
    raise Exception(f"Batch transform job {error} failed.")

print("Next step")


# # Smoothing Layer Inference


jobs = []
for SUBMODEL_NAME in SUBMODELS:
    print(f"Running submodel {SUBMODEL_NAME}")

    job_name = f"smooth-{uuid.uuid4()}"
    jobs.append(job_name)

    # get directory
    SMOOTH_OUTPUT_DIR = f"{OUTPUT_DIR.smooth_results}{SUBMODEL_NAME}/"
    print(SMOOTH_OUTPUT_DIR)

    # create model
    model_registry = get_model_registry(ANCESTRY_MODEL, [SUBMODEL_NAME])

    model = PyTorchModel(
        model_data=model_registry["sub-models"][SUBMODEL_NAME]["model"][
            "smooth_model_uri"
        ],
        role=ROLE,
        framework_version="1.10",
        py_version="py38",
        source_dir="src/",
        entry_point="inference_user.py",
    )

    # batch transformer
    transformer = model.transformer(
        output_path=SMOOTH_OUTPUT_DIR,
        instance_count=1,
        instance_type="ml.m5.large",
        accept="text/tsv",
        strategy="MultiRecord",
        max_payload=4,
    )

    # run transformer
    BASE_OUTPUT_DIR = f"{OUTPUT_DIR.base_results}{SUBMODEL_NAME}/"
    transformer.transform(
        data=BASE_OUTPUT_DIR,
        content_type="text/tsv",
        split_type="Line",
        job_name=job_name,
        model_client_config={"InvocationsMaxRetries": 0},
        logs=False,
        wait=False,
    )


client = boto3.client("sagemaker")
errors = []

while len(jobs):
    print(f"{len(jobs)} left")
    for job_name in jobs:
        response = client.describe_transform_job(
            TransformJobName=job_name,
        )["TransformJobStatus"]

        if response in ["InProgress", "Stopping"]:
            time.sleep(60)
            continue

        if response in ["Failed", "Stopped"]:
            errors.append(job_name)

        indx = jobs.index(job_name)
        jobs.pop(indx)
        break

for error in errors:
    raise Exception(f"Smooth transform job {error} failed.")

print("Next step")


# # Base & smoothing layers post-processing


# convert base layer results to summary format
# generate summary of the smoothing layer output
for SUBMODEL_NAME in SUBMODELS:
    print(f"Running submodel {SUBMODEL_NAME}")

    # model registry
    model_registry = get_model_registry(ANCESTRY_MODEL, [SUBMODEL_NAME])

    parameters_file = model_registry["sub-models"][SUBMODEL_NAME]["artifacts"][
        "parameters"
    ]

    os.system("rm -rf temp/")
    os.system("mkdir -p temp/")
    os.system(f"aws s3 cp {parameters_file} temp/ --quiet")

    with open(f"temp/parameters.json", "rt") as fin:
        params = json.load(fin)

    for layer in ["base", "smooth"]:
        # base output path
        path = getattr(OUTPUT_DIR, f"{layer}_results")
        LAYER_OUTPUT_DIR = f"{path}{SUBMODEL_NAME}/"

        # read data and convert to regular format
        os.system(
            f"aws s3 cp --recursive {LAYER_OUTPUT_DIR} temp/{layer}/intermediate/ --quiet"
        )

        # read base results & save
        df_base = convert_to_ancestry_format(
            dd.read_csv(
                f"temp/{layer}/intermediate/*",
                sep="\t",
                header=None,
                dtype={2: "object"},
            ).compute(),
            model_registry["artifacts"]["population_map_uri"],
            params,
            pred_by_argmin=bool(layer == "base"),
        )
        df_base.to_csv(
            f"{OUTPUT_DIR.base}{OUTPUT_DIR.path}summary-results/{layer}_samples.{SUBMODEL_NAME}.tsv.gz",
            sep="\t",
            index=False,
        )
        del df_base


# ### Performing basic analysis of the results (saved to ancestry.tsv)


os.system("mkdir -p temp/")
os.system(
    f"aws s3 cp {OUTPUT_DIR.base}{OUTPUT_DIR.path}summary-results/ temp/ --recursive"
)

df_base = dd.read_csv(
    "temp/base_samples.*.tsv.gz",
    sep="\t",
    blocksize=None,
    dtype={"sample_id": "object"},
).compute()
df_base["layer"] = "base"

df_smooth = dd.read_csv(
    "temp/smooth_samples.*.tsv.gz",
    sep="\t",
    blocksize=None,
    dtype={"sample_id": "object"},
).compute()
df_smooth["layer"] = "smooth"

df = pd.concat([df_smooth, df_base])
del df_base, df_smooth

with open("chr_map.json", "rt") as fin:
    chr_fractions = json.load(fin)

chr_fractions = {int(k): v for k, v in chr_fractions.items()}
window_lengths = df.groupby(["chrom"])["window"].nunique().to_dict()
chr_map = {
    k: v / window_lengths[k] for k, v in chr_fractions.items() if k in window_lengths
}

df["weight"] = df.chrom.map(chr_map)


def weight_to_str(x):
    x = x / x.sum()
    x = x.sort_values(ascending=False)
    return ", ".join(f"{idx[2]}: {v*100:0.2f}" for idx, v in x.items())


df_ = (
    df.groupby(["sample_id", "layer", "pred"])
    .agg({"weight": "sum"})
    .groupby(["sample_id", "layer"])
    .agg({"weight": weight_to_str})
    .rename(columns={"weight": "predicted ancestry"})
    .reset_index()
    .sort_values(["sample_id", "layer"])
)
df_.to_csv(f"{NAME}-ancestry.tsv", sep="\t", index=False)

os.system("rm -rf temp/")
