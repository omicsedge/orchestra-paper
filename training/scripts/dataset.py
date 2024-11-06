import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset


class LaiDataset(Dataset):
    """LAI Dataset for Smoothing Layer."""

    def __init__(self, path: str, device) -> None:
        super().__init__()

        self._path = Path(path)
        self._device = device

        _, X, self.n_populations = self.read_X(self._path)
        self._X = torch.from_numpy(X)

        y = self.read_y(self._path)
        self._y = torch.from_numpy(y)

    @classmethod
    def read_X(cls, path: str):
        df = pd.read_csv(path / "predictions.tsv.gz", sep="\t")

        n_populations, windows = cls._get_parameters(path)
        pred_cols = cls._get_pred_columns(n_populations, windows)

        X = df.iloc[:, 1 : 1 + pred_cols].to_numpy().astype(np.float32)
        n_samples = X.shape[0]
        n_windows = pred_cols // n_populations
        X = X.reshape(n_samples, n_windows, n_populations)

        sample_ids = df.sample_id.to_numpy()

        return sample_ids, np.swapaxes(X, 1, 2), n_populations

    @classmethod
    def read_y(cls, path: str):
        df = pd.read_csv(path / "predictions.tsv.gz", sep="\t")

        n_populations, windows = cls._get_parameters(path)
        pred_cols = cls._get_pred_columns(n_populations, windows)

        return df.iloc[:, 1 + pred_cols :].to_numpy().astype(np.int64)

    @classmethod
    def _get_pred_columns(self, n_populations, windows):
        pred_cols = 0
        for w in windows.values():
            pred_cols += w * n_populations

        return pred_cols

    @classmethod
    def _get_parameters(cls, path):
        with open(path / "parameters.json", mode="rt") as fin:
            params = json.load(fin)
        n_populations = params["n-populations"]
        windows = params["windows"]  # windows per chromosome
        return n_populations, windows

    def __len__(self):
        return len(self._y)

    def __getitem__(self, idx):
        return {
            "inputs": self._X[idx].to(self._device),
            "labels": self._y[idx].to(self._device),
        }

    @property
    def X(self):
        return self._X

    @property
    def y(self):
        return self._y
