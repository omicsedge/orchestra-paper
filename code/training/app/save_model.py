import logging
import os
import shutil
import tarfile
from pathlib import Path

logging.getLogger().setLevel(logging.INFO)


class AncestryModel:
    def __init__(
        self,
        population_map_file: str,
        model_path: Path,
        sub_model_name: str,
        windows_info_file: Path,
        base_model_path: Path,
        smooth_model_path: str,
    ) -> None:
        self.model_path = model_path
        self.population_map_file = population_map_file
        self.windows_info_file = windows_info_file
        self.base_model_path = base_model_path
        self.smooth_model_path = smooth_model_path
        self.submodel_path = model_path / "sub-models" / sub_model_name

    def _create_dirs(self):
        (self.model_path / "artifacts").mkdir(parents=True, exist_ok=True)
        (self.submodel_path / "artifacts").mkdir(parents=True, exist_ok=True)
        (self.submodel_path / "base").mkdir(parents=True, exist_ok=True)
        (self.submodel_path / "smooth").mkdir(parents=True, exist_ok=True)

    def _save_artifacts(self):
        shutil.copyfile(
            self.population_map_file,
            self.model_path / "artifacts" / "population_map.tsv",
        )
        shutil.copyfile(
            self.windows_info_file,
            self.submodel_path / "artifacts" / "windows_info.tsv",
        )
        shutil.copyfile(
            self.smooth_model_path / "parameters.json",
            self.submodel_path / "artifacts" / "parameters.json",
        )

    def _make_tarfile(self, source_dir: Path, output_filename: Path):
        # Ensure the source directory exists
        if not os.path.exists(source_dir):
            logging.error(f"Directory '{source_dir}' does not exist!")
            return

        with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname="model")

    def _save_sub_model(self):
        self._make_tarfile(
            self.base_model_path, self.submodel_path / "base" / "model.tar.gz"
        )
        self._make_tarfile(
            self.smooth_model_path, self.submodel_path / "smooth" / "model.tar.gz"
        )
        logging.info("Models created successfully")

    def save(self):
        self._create_dirs()
        self._save_artifacts()
        self._save_sub_model()
