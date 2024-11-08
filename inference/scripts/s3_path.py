from urllib.parse import urlparse
import os


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
