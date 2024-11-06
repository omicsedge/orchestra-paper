import logging
import os
import shutil
import tarfile
from pathlib import Path

logging.getLogger().setLevel(logging.INFO)


class AncestryModel:
    def __init__(
        self,
        simulated_data_path: Path,
        model_path: Path,
        sub_model_name: str,
        windows_info_file: Path,
        simulated_params_file: str,
        base_model_path: Path,
        smooth_model_path: str,
    ) -> None:
        self.model_path = model_path
        self.simulated_data_path = simulated_data_path
        self.windows_info_file = windows_info_file
        self.simulated_params_file = simulated_params_file
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
            self.simulated_data_path / "population_map.tsv",
            self.model_path / "artifacts" / "population_map.tsv",
        )
        shutil.copyfile(
            self.windows_info_file,
            self.submodel_path / "artifacts" / "windows_info.tsv",
        )
        shutil.copyfile(
            self.simulated_params_file,
            self.submodel_path / "artifacts" / "parameters.json",
        )

    def _make_tarfile(self, output_filename, source_dir):
        # Ensure the source directory exists
        if not os.path.exists(source_dir):
            logging.error(f"Directory '{source_dir}' does not exist!")
            return

        with tarfile.open(output_filename, "w:gz") as tar:
            tar.add(source_dir, arcname=os.path.basename(source_dir))

    def _save_sub_model(self):
        self._make_tarfile(
            self.submodel_path / "base" / "model.tar.gz", self.base_model_path
        )
        shutil.copyfile(
            self.smooth_model_path,
            self.submodel_path / "smooth" / "model.pt",
        )
        logging.info("Models created successfully")

    def save(self):
        self._create_dirs()
        self._save_artifacts()
        self._save_sub_model()
