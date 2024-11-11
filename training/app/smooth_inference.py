import json
import logging
import os
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import torch

logging.basicConfig(level=logging.DEBUG)


def read_x(path: Path, n_populations: int, windows: Dict[str, int]):
    df = pd.read_csv(path, sep="\t")
    pred_cols = _get_pred_columns(n_populations, windows)

    X = df.iloc[:, 1 : 1 + pred_cols].to_numpy().astype(np.float32)
    n_samples = X.shape[0]
    n_windows = pred_cols // n_populations
    X = X.reshape(n_samples, n_windows, n_populations)

    sample_ids = df.sample_id.to_numpy()

    return sample_ids, np.swapaxes(X, 1, 2)


def read_y(path: Path, n_populations: int, windows: Dict[str, int]):
    df = pd.read_csv(path, sep="\t")
    pred_cols = _get_pred_columns(n_populations, windows)
    return df.iloc[:, 1 + pred_cols :].to_numpy().astype(np.int64)


def _get_pred_columns(n_populations, windows):
    pred_cols = 0
    for w in windows.values():
        pred_cols += w * n_populations

    return pred_cols


def model_fn(model_dir: Path):
    model = torch.jit.load(model_dir / "model.pt", map_location=torch.device("cpu"))
    model.eval()
    return model


def predict_fn(data, model):
    batch_size = 1024
    device = torch.device("cpu")

    X = torch.from_numpy(data).to(device)
    y = np.zeros_like(data)

    n_batches = (
        X.shape[0] // batch_size
        if X.shape[0] % batch_size == 0
        else X.shape[0] // batch_size + 1
    )
    with torch.no_grad():
        for i in range(n_batches):
            logging.info(f"Batch {i+1}")
            X_ = X[i * batch_size : (i + 1) * batch_size, :, :]
            y_pred = model(X_)
            y_pred = torch.softmax(y_pred, dim=1)
            y[i * batch_size : (i + 1) * batch_size, :, :] = y_pred.numpy()

    return np.swapaxes(y, 1, 2)


def sample_id_to_gen(sid):
    if "generation" not in sid:
        return 0
    else:
        return int(sid.split("_")[-2][len("generation") :])


def save_results(sample_ids, y_pred, y_true, windows: Dict[str, int], output_dir: Path):
    n_pops = y_pred.shape[2]

    generations = np.vectorize(sample_id_to_gen)(sample_ids)

    start_window = 0
    for chrom, n_windows in sorted(windows.items(), key=lambda x: x[0]):
        print(f"Saving chromosome {chrom}")

        predict_dir = output_dir / f"chr{chrom}"
        os.makedirs(predict_dir, exist_ok=True)

        for i in range(n_windows):
            w = i + start_window

            for g in np.unique(generations):
                idx = generations == g

                df = pd.concat(
                    [
                        pd.DataFrame(sample_ids[idx], columns=["sample_id"]),
                        pd.DataFrame(
                            y_pred[idx, w, :], columns=[f"{p}" for p in range(n_pops)]
                        ),
                        pd.DataFrame(y_true[idx, w], columns=["label"]),
                    ],
                    axis=1,
                )

                output_file = (
                    predict_dir / f"gen{g}" / f"window-{i}" / "prediction.tsv.gz"
                )
                os.makedirs(output_file.parents[0], exist_ok=True)
                df.to_csv(output_file, sep="\t", index=False)
                logging.info(f"saved {output_file} (n_samples: {len(df)}")

        start_window = start_window + n_windows


def smooth_inference(
    smooth_layer_model: Path,
    base_layer_pred: Path,
    base_layer_params: Path,
    output_dir: Path,
):
    print("start running....")
    output_dir.mkdir(parents=True, exist_ok=True)

    # output the results
    with open(base_layer_params, mode="rt") as fin:
        params = json.load(fin)
        n_populations = params["n-populations"]
        windows = params["windows"]

    # load input data
    sample_ids, X = read_x(base_layer_pred, n_populations, windows)
    y_true = read_y(base_layer_pred, n_populations, windows)

    logging.info(f"Number of samples: {len(sample_ids)}")
    logging.info(f"Number of windows: {X.shape[1]}")
    logging.info(f"Number of populations: {X.shape[2]}")

    # predict
    model = model_fn(model_dir=smooth_layer_model)
    y_pred = predict_fn(X, model)

    save_results(sample_ids, y_pred, y_true, windows, output_dir)
