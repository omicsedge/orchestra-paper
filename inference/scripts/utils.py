import logging

import numpy as np
import pandas as pd

# configure logging
logging.basicConfig(level=logging.INFO)


def convert_to_ancestry_format(
    df, population_map_file, parameters, pred_by_argmin=False
):
    df_map = pd.read_csv(population_map_file, sep="\t")
    id_to_code = df_map.set_index("id")["level_3_code"].to_dict()

    n_populations = len(df_map)  # optional len(set(id_to_code.values()))
    logging.info(f"Num populations: {n_populations}")

    n_haps = len(df)
    logging.info(f"Num haplotypes: {n_haps}")

    n_windows = sum(parameters["windows"].values())
    logging.info(f"Number of windows: {n_windows}")

    y = df.iloc[:, 4:].to_numpy()
    if n_populations * n_windows != y.shape[1]:
        raise ValueError(
            f"Mismatch in number of columns. "
            f"Expected {n_populations * n_windows}, but received {y.shape[1]}"
        )

    # reshape predictions into per window for raw format
    y = y.reshape(y.shape[0] * n_windows, n_populations)
    logging.info(f"Shape of raw predictiosn: {y.shape}")

    # create dataframe
    df_hap = pd.DataFrame(y)
    df_hap.columns = (
        df_map.level_3_code.tolist()
    )  # optional list(set(df_map.level_3_code.tolist()))

    # add file id
    file_id = np.repeat(df.iloc[:, 0].tolist(), n_windows)
    df_hap.insert(0, "file_id", file_id)

    # add sample_id
    sample_id = np.repeat(df.iloc[:, 2].tolist(), n_windows)
    df_hap.insert(1, "sample_id", sample_id)

    # add haplotype index
    hap = np.repeat(df.iloc[:, 3].tolist(), n_windows)
    df_hap.insert(2, "hap", hap)

    # add chromosome x windows
    chroms = []
    windows = []
    for chrom, nw in parameters["windows"].items():
        chroms += [int(chrom)] * nw
        windows += list(range(nw))
    df_hap.insert(3, "chrom", np.tile(chroms, n_haps))
    df_hap.insert(4, "window", np.tile(windows, n_haps))

    # add prediction
    pred_fn = np.argmin if pred_by_argmin else np.argmax
    df_hap.insert(5, "pred", pred_fn(y, axis=1))
    df_hap["pred"] = df_hap.pred.map(id_to_code)

    return df_hap
