import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import softmax
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score

logging.basicConfig(level=logging.DEBUG)


def get_data_from_path(path: Path, data_name: str):
    return path.parts[-1][len(data_name) :]


def save_evaluation(filename, metrics):
    with open(filename, mode="wt") as fout:
        json.dump(metrics, fout, indent=4)


def calculate_metrics(df):
    n_populations = len([c for c in df.columns if c.isdigit()])
    logging.info(f"Number of Populations: {n_populations}")

    # get true labels
    y_true = df.label.to_numpy()

    # get predictions
    preds = df.iloc[:, 1 : 1 + n_populations].to_numpy()
    y_pred = np.argmax(preds, axis=1)
    df["pred"] = y_pred

    # calculate accuracy
    accuracy = accuracy_score(y_true, y_pred)
    logging.info(f"Accuracy: {accuracy*100:0.2f}%")

    # calculate accuracy by generation - fixed column selection syntax
    gen_accuracy = (
        df.groupby("gen")[["label", "pred"]]
        .apply(lambda x: np.mean(x.label == x.pred))
        .to_dict()
    )

    # generate output dictionary
    output = {
        "multiclass_classification_metrics": {
            "accuracy": {
                "value": accuracy,
            },
            "weighted_recall": {
                "value": recall_score(y_true, y_pred, average="weighted"),
            },
            "weighted_precision": {
                "value": precision_score(y_true, y_pred, average="weighted"),
            },
            "weighted_f1": {
                "value": f1_score(y_true, y_pred, average="weighted"),
            },
            # "confusion_matrix": cm_dict,
            "generations": gen_accuracy,
        },
    }

    return output


def load_recomb_files(input_dir: Path, is_recomb):
    files = input_dir.rglob("prediction.tsv.gz")

    dfs = []
    for _file in files:
        df = pd.read_csv(_file, sep="\t")
        if is_recomb:
            # Apply softmax with more stable computation
            # df.iloc[:, 1:-1] = softmax(1 / np.clip(df.iloc[:, 1:-1], a_min=1e-10, a_max=None), axis=1)
            df.iloc[:, 1:-1] = softmax(1 / df.iloc[:, 1:-1], axis=1)
        df["window"] = int(_file.parent.name.split("-")[-1])
        df["gen"] = int(_file.parent.parent.name[len("gen") :])
        dfs.append(df)

    df = pd.concat(dfs)
    return df


def run(input_dir: Path, output_dir: Path, is_recomb: bool):
    logging.info(f"Input Dir: {input_dir}")

    chrom = int(get_data_from_path(input_dir, "chr"))
    logging.info(f"Chromosome: {chrom}")

    df = load_recomb_files(input_dir, is_recomb)

    n_samples = df.sample_id.nunique()
    n_windows = df.window.nunique()
    logging.info(f"Number of Samples: {n_samples}")
    logging.info(f"Number of Windows: {n_windows}")

    output = calculate_metrics(df)
    logging.info(json.dumps(output, indent=4))

    output_dir.mkdir(parents=True, exist_ok=True)
    save_evaluation(output_dir / f"chr{chrom}.json", output)

    return df


def evaluate(input_dir: Path, output_dir: Path, is_recomb: bool = False):
    output_dir.mkdir(parents=True, exist_ok=True)

    logging.info(f"input directory: {input_dir}")
    logging.info(f"output directory: {output_dir}")
    dfs = []

    for chrom_dir in input_dir.glob("chr*"):
        df = run(chrom_dir, output_dir, is_recomb)
        dfs.append(df)

    output = calculate_metrics(pd.concat(dfs))
    logging.info(json.dumps(output, indent=4))

    save_evaluation(output_dir / "all.json", output)
