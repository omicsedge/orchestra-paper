# get stitched samples (recomb distance) and its dependencies uses only numpy
import json
import logging
import os
from io import BytesIO, StringIO
from pathlib import Path
from typing import Tuple

import dask.dataframe as dd
import numpy as np
import pandas as pd
import torch
from dask import delayed as ddelayed
from joblib import Parallel, delayed
from numba import njit
from scipy.special import softmax

# configure logging
logging.basicConfig(level=logging.INFO)


BATCH_SIZE = 1024

class Recombination:
    def __init__(self):
        pass

    def fit(self, X, y):
        self.n_populations = len(np.unique(y))
        self.pops = np.arange(self.n_populations)
        self.X_reference = X
        self.y_reference = y

        return self

    def transform(self, X, train_set):
        # store the n_segments of each sample against each population
        recombination_dist_output = np.empty((X.shape[0], self.n_populations))

        # iterate over samples in input
        for sample_id, sample_data in enumerate(X):
            # get similarity matrix of full reference population
            similarity_matrix = self.X_reference == sample_data
            y = self.y_reference

            if train_set:
                # remove parents
                similarity_matrix = remove_parents(similarity_matrix)

            # calculate n_segments for each population
            for pop in self.pops:

                # select only samples in reference  reference population
                select_population = (y == pop).ravel()

                # get recombination distance
                recombination_dist_output[sample_id][pop] = get_recombination_distance(
                    similarity_matrix[select_population]
                )

        return recombination_dist_output

    def predict_proba(self, X, train_set, dimension=1):
        base_pred = self.transform(X, train_set)
        normalized = (
            np.expand_dims(np.apply_along_axis(np.max, dimension, base_pred), dimension)
            - base_pred
        )
        return np.apply_along_axis(softmax, dimension, normalized)


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


class Transformer:
    """Transform imput using recomb base layer."""

    def __init__(self, model_dir: Path):
        """Initialize the transformer."""
        self._model_dir = model_dir

    def transform(self, data: pd.DataFrame, output_file: str, n_jobs=None):
        """
        Make predictions for input data using recomb base layer.

        Args:
            data pd.DataFrame: genotype data for the samples
            output_file str: output file path
            n_jobs (_type_, optional): number of parallel jobs to run

        Returns:
            numpy array with shape `n_samples x (n_populations x n_windows)`
        """

        if n_jobs is None:
            n_jobs = os.cpu_count()

        # parse data
        df_labels = data.iloc[:, :4]
        genotypes = data.iloc[:, 4:].to_numpy(dtype=np.int64)

        windows = self._load_model_windows()
        window_data = np.hsplit(genotypes, np.cumsum(windows["length"].to_numpy()))

        res = Parallel(n_jobs=n_jobs)(
            delayed(self.transform_window)(chrom, w, X)
            for (chrom, w), X in zip(
                windows[["chr", "window"]].itertuples(index=False, name=None),
                window_data,
            )
        )

        predictions = np.hstack(res)
        print("Prediction shape: ", predictions.shape)

        # prepare response
        df = pd.concat([df_labels, pd.DataFrame(predictions)], axis=1)
        df.to_csv(output_file, header=False, index=False, sep="\t")

    def transform_window(self, chrom: int, window: int, data: np.ndarray):
        X_reference, y_reference = self._load_model(chrom, window)

        model = Recombination()
        model.fit(X_reference, y_reference)

        return model.transform(data, train_set=False)

    def _load_model(self, chrom: int, window: int) -> Tuple[np.ndarray, np.ndarray]:
        model_file = self._model_dir / f"chr{chrom}/window-{window}/data.parquet"
        df_model = pd.read_parquet(model_file)
        x_reference = df_model.iloc[:, 3:-1].to_numpy()
        y_reference = df_model["label"].to_numpy()

        return x_reference, y_reference

    def _load_model_windows(self) -> pd.DataFrame:
        snp_files = self._model_dir.rglob("chr*/window-*/snps.tsv")
        dfs = [ddelayed(self._read_window(snp_file)) for snp_file in snp_files]
        df = dd.from_delayed(dfs).compute()
        return df.groupby(["chr", "window"])["pos"].count().reset_index(name="length")

    @classmethod
    def _read_window(cls, filename) -> pd.DataFrame:
        # reads each csv file to a pandas.DataFrame
        df = pd.read_csv(filename, sep="\t")
        df["window"] = int(filename.parts[-2][len("window-") :])
        return df

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

    if "text/tsv" in request_content_type:
        df = pd.read_csv(
            StringIO(BytesIO(request_body).read().decode()), sep="\t", header=None
        )
    else:
        logging.error(f"not support {request_content_type} type yet")
        raise ValueError(f"not support {request_content_type} type yet")

    return df




@njit
def get_recombination_distance(similarity_matrix):
    """
    Parameters
    ----------
    similarity_matrix : boolean 2D numpy array (n_reference_samples,window_size)
        DESCRIPTION.
    Returns
    -------
    TYPE int
        recombination distance.
    """
    _, nc = similarity_matrix.shape

    rd = 0
    solution = []
    index_tracker = -1
    pointer = 0
    start = 0

    running_product = True + similarity_matrix[:, 0]

    for pointer in range(nc):
        running_product *= similarity_matrix[:, pointer]

        if running_product.max() == 0:
            rd += 1
            solution.append((index_tracker, start, pointer))
            start = pointer
            # running_product =  True + similarity_matrix[:,0] #no overlap, cause we want to remove parents.
            running_product = similarity_matrix[:, pointer]

        index_tracker = np.argmax(running_product)

    rd += 1
    solution.append((index_tracker, start, pointer))

    for ind, start, stop in solution:
        similarity_matrix[ind, start:stop] = False

    rd1 = rd

    rd = 0
    running_product = True + similarity_matrix[:, 0]

    for pointer in range(nc):
        running_product *= similarity_matrix[:, pointer]

        if running_product.max() == 0:
            rd += 1
            running_product = similarity_matrix[:, pointer]  # enforce overlap of one

    return rd1 * rd / (rd1 + rd)


@njit
def remove_parents(similarity_matrix):
    _, nc = similarity_matrix.shape

    rd = 0
    solution = []
    index_tracker = -1
    pointer = 0
    start = 0

    running_product = True + similarity_matrix[:, 0]

    for pointer in range(nc):
        running_product *= similarity_matrix[:, pointer]

        if running_product.max() == 0:
            rd += 1
            solution.append((index_tracker, start, pointer))
            start = pointer
            running_product = (
                True + similarity_matrix[:, 0]
            )  # no overlap, cause we want to remove parents.

        index_tracker = np.argmax(running_product)

    rd += 1
    solution.append((index_tracker, start, pointer))
    # solution = np.array(solution, dtype = int)

    for ind, start, stop in solution:
        similarity_matrix[ind, start:stop] = False

    return similarity_matrix


def predict_fn(df, model_path):
    logging.debug("**** Model Function *****")
    model = Model(model_path)
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
