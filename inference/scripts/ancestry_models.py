import json
import logging
import os
import tarfile
import tempfile
import typing
from urllib.parse import urlparse

logger = logging.getLogger()
logger.setLevel(logging.INFO)


def s3_uri_to_parts(s3_uri):
    o = urlparse(s3_uri)
    return o.netloc, o.path.lstrip("/")


def s3_copy(from_uri, to_uri):
    logging.info(f"copy {from_uri} to {to_uri}")
    bucket, key = s3_uri_to_parts(from_uri)
    copy_source = {
        "Bucket": bucket,
        "Key": key,
    }

    s3 = boto3.resource("s3")

    bucket, key = s3_uri_to_parts(to_uri)
    s3.meta.client.copy(copy_source, Bucket=bucket, Key=key)


def s3_download(from_uri, to_file):
    logging.info(f"download {from_uri} to {to_file}")
    bucket, key = s3_uri_to_parts(from_uri)

    s3 = boto3.client("s3")
    s3.download_file(bucket, key, to_file)


def s3_upload_fileobj(fin, to_uri):
    bucket, key = s3_uri_to_parts(to_uri)

    s3 = boto3.client("s3")
    s3.upload_fileobj(fin, bucket, key)


class AncestrySubModel:
    def __init__(self, model: "AncestryModel", name: str) -> None:
        self._model = model
        self._name = name

    @property
    def model_path(self):
        return self._model.join([self._model.model_path, "sub-models", self._name])

    @property
    def base_model_uri(self):
        return self._model.join([self.model_path, "base/model.tar.gz"])

    @property
    def smooth_model_uri(self):
        return self._model.join([self.model_path, "smooth/model.tar.gz"])

    @property
    def parameter_file_uri(self):
        return self._model.join([self.model_path, "artifacts/parameters.json"])

    @property
    def window_info_file_uri(self):
        return self._model.join([self.model_path, "artifacts/windows_info.tsv"])

    def save(
        self,
        base_model_uri,
        base_output_uri,
        smooth_model_uri,
        parameter_path,
        dataset_path,
    ):
        s3_copy(base_model_uri, self.base_model_uri)
        s3_copy(smooth_model_uri, self.smooth_model_uri)
        s3_copy(
            os.path.join(parameter_path, "parameters.json"), self.parameter_file_uri
        )
        s3_copy(
            os.path.join(dataset_path, "population_map.tsv"),
            self._model.population_map_uri,
        )

        with tempfile.TemporaryDirectory() as tmpdirname:
            output_file = os.path.join(tmpdirname, "output.tar.gz")
            s3_download(base_output_uri, output_file)
            tar = tarfile.open(output_file, mode="r:gz")
            for tar_info in tar:
                if tar_info.name != "windows_info.tsv":
                    continue

                file_save = tar.extractfile(tar_info.name)
                s3_upload_fileobj(file_save, self.window_info_file_uri)

            tar.close()

    @property
    def chromosomes(self) -> str:
        """Get chromosome list for sql query"""
        return "'" + "', 'chr".join(self._name.split(".")) + "'"


class AncestryModel:
    def __init__(
        self,
        model_registry_uri,
        simulated_dataset,
        version,
        is_pipeline=False,
    ) -> None:
        self._is_pipeline = is_pipeline
        self._model_registry_uri = model_registry_uri
        self._simulated_dataset = simulated_dataset
        self._version = version
        self._submodels = self._get_submodels()

    def join(self, path_segments: typing.List[str]):
        return os.path.join(*path_segments)

    @property
    def name(self):
        return self.join([self._simulated_dataset, self._version])

    @property
    def model_path(self):
        return self.join(
            [self._model_registry_uri, self._simulated_dataset, self._version]
        )

    @property
    def population_map_uri(self):
        return self.join([self.model_path, "artifacts/population_map.tsv"])

    @property
    def submodels(self):
        return self._submodels

    def get_submodel(self, name):
        return AncestrySubModel(model=self, name=name)

    def _get_submodels(self):
        bucket, key = s3_uri_to_parts(self.model_path)
        client = boto3.client("s3")
        result = client.list_objects_v2(
            Bucket=bucket, Prefix=key + "/sub-models/", Delimiter="/"
        )
        return [o.get("Prefix").split("/")[-2] for o in result.get("CommonPrefixes")]


def save_model_lambda_handler(event, context):
    logging.info(json.dumps(event, indent=4))
    model = AncestryModel(
        event["model-registry-uri"], event["simulated-dataset"], event["version"]
    )
    sub_model = model.get_submodel(event["sub-model-name"])
    sub_model.save(
        base_model_uri=event["base-model-uri"],
        base_output_uri=event["base-output-uri"],
        smooth_model_uri=event["smooth-model-uri"],
        parameter_path=event["parameter-file-path"],
        dataset_path=event["dataset-path"],
    )

    return {
        "statusCode": 200,
    }
