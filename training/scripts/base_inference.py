"""
get stitched samples (recomb distance) and its dependencies uses only numpy
"""

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit
from scipy.special import softmax

# configure logging
logging.basicConfig(level=logging.INFO)


def start_logger():
    logger = logging.getLogger("mylogger")
    if len(logger.handlers) == 0:
        logger.setLevel(logging.INFO)
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(message)s"))
        logger.addHandler(sh)
    return logger


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


def get_item_from_path(path: Path, prefix: str):
    item = None
    for part in path.parts:
        if part.startswith(prefix):
            item = int(part[len(prefix) :])
            break

    if item is None:
        raise FileNotFoundError()

    return item


def predict_window(
    data_file: Path,
    base_layer_model: Path,
    output_dir: Path,
    training_set: bool,
):
    gen = get_item_from_path(data_file, "gen")
    window = get_item_from_path(data_file, "window-")
    chrom = get_item_from_path(data_file, "chr")

    logger = start_logger()
    logger.info(f"predicting window chr{chrom}:gen{gen}:{window} ({data_file})...")

    # load reference panel
    reference_file = base_layer_model / f"window-{window}" / "data.parquet"

    df_reference = pd.read_parquet(reference_file)
    x_reference = df_reference.iloc[:, 3:-1].to_numpy()
    y_reference = df_reference["label"].to_numpy()

    # load test data
    df = pd.read_parquet(data_file)
    X = df.iloc[:, 3:-1].to_numpy()

    # predict
    model = Recombination()
    model.fit(x_reference, y_reference)
    y_pred = model.transform(X, train_set=training_set)

    df_pred = pd.concat(
        [
            pd.DataFrame({"sample_id": df.sample_id}),
            pd.DataFrame(y_pred),
        ],
        axis=1,
    )

    # calculate accuracy for
    if "label" in df.columns:
        y = df["label"].to_numpy()
        df_pred["label"] = y

        y_ = np.argmin(y_pred, 1)
        accuracy = (y_ == y.ravel()).mean()
        logger.info(
            f"Window chr{chrom}:gen{gen}:{window:3} accuracy: {accuracy * 100:0.2f}%"
        )

    # save output
    output_path = output_dir / f"gen{gen}" / f"window-{window}"
    output_path.mkdir(parents=True, exist_ok=True)
    df_pred.to_csv(output_path / "prediction.tsv.gz", sep="\t", index=False)
    logger.info(f"Output path: {output_path}")


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
