# get stitched samples (recomb distance) and its dependencies uses only numpy
import itertools
import json
import logging
import os
from importlib.resources import path
from io import BytesIO, StringIO
from pathlib import Path

import click
import numpy as np
import pandas as pd
import torch
from sagemaker_inference import encoder

logging.basicConfig(level=logging.DEBUG)


BATCH_SIZE = 1024


class Model:
    def __init__(self, model_dir):
        model_file = os.path.join(model_dir, "model.pt")
        self._model = torch.jit.load(model_file, map_location=torch.device("cpu"))
        self._model.eval()

        parameter_file = os.path.join(model_dir, "parameters.json")
        with open(parameter_file) as fin:
            self._parameters = json.load(fin)

    @property
    def model(self):
        return self._model

    @property
    def n_populations(self):
        return self._parameters["n-populations"]

    @property
    def windows_per_chromosome(self):
        return self._parameters["windows"]

    def predict_proba(self, X):
        return self._model(X)


def model_fn(model_dir):
    logging.debug("**** Model Function *****")
    return Model(model_dir)


def read_X(df, n_populations):
    X = df.iloc[:, 4:].to_numpy().astype(np.float32)
    n_samples = X.shape[0]
    n_windows = X.shape[1] // n_populations
    logging.info(
        f"n_samples x n_windows x n_populations: {n_samples} x {n_windows} x {n_populations}"
    )

    assert X.shape[1] == n_populations * n_windows
    X = X.reshape(n_samples, n_windows, n_populations)

    return np.swapaxes(
        X, 1, 2
    )  # convert to (n_samples, n_populations, n_windows) shape


def input_fn(request_body, request_content_type):
    logging.debug("**** Input Function *****")

    # TODO: check whether text/tsv is in the content type
    if "text/tsv" in request_content_type:
        df = pd.read_csv(
            StringIO(BytesIO(request_body).read().decode()), sep="\t", header=None
        )
    else:
        logging.error(f"not support {request_content_type} type yet")
        raise ValueError(f"not support {request_content_type} type yet")

    return df


def predict_fn(df, model):
    logging.debug("**** Predict Function *****")
    try:
        device = torch.device("cpu")

        X = read_X(df, model.n_populations)
        y = np.zeros_like(X)

        logging.debug("X shape: ", X.shape)
        X = torch.from_numpy(X).to(device)

        n_batches = (
            X.shape[0] // BATCH_SIZE
            if X.shape[0] % BATCH_SIZE == 0
            else X.shape[0] // BATCH_SIZE + 1
        )
        logging.debug(f"Number of batches: {n_batches}")

        with torch.no_grad():
            for i in range(n_batches):
                logging.info(f"Batch {i+1}/{n_batches}")
                X_ = X[i * BATCH_SIZE : (i + 1) * BATCH_SIZE, :, :]
                logging.info(X_.shape)

                logging.info("  predicting...")
                y_pred = model.predict_proba(X_)
                y_pred = torch.softmax(y_pred, dim=1)
                logging.info(f"  predicted: shape({y_pred.shape})...")

                y[i * BATCH_SIZE : (i + 1) * BATCH_SIZE, :, :] = y_pred.numpy()
                logging.info("  batch added.")

        logging.info(f"prepare output data: {y.shape}")
        y = np.swapaxes(y, 1, 2).reshape(y.shape[0], -1)

        df_output = pd.DataFrame(y)
        df_output = pd.concat([df.iloc[:, :4], df_output], axis=1)

        logging.info(df_output.head())
    except Exception as exp:
        logging.error(exp)
        raise exp

    return df_output


def output_fn(df, accept):
    """A default output_fn for PyTorch. Serializes predictions from predict_fn to JSON, CSV or NPY format.

    Args:
        df: a prediction result from predict_fn
        accept: type which the output data needs to be serialized

    Returns: output data serialized
    """
    logging.info(f"accept: {accept}")
    out = StringIO()
    df.to_csv(out, header=False, index=False, sep="\t")

    # return encoder.encode(out.getvalue(), accept)
    return out.getvalue()
