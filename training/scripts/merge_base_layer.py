import csv
import gzip
import json
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.DEBUG)


class GenerationDataFrameIterator:
    def __init__(self, input_dir: Path):
        self._window_info = {}
        self._gen_df_itrs = self._get_generation_dataframe_iterators(input_dir)
        self._n_populations = None

        self._gen_dfs = {}
        self._load_dfs()
        self._finished = False

    @property
    def window_info(self):
        return self._window_info

    @property
    def n_populations(self):
        return self._n_populations

    def get_next_batch(self):
        if self._finished:
            return None

        rows = []
        should_stop = False
        i = 0
        while not should_stop:
            should_stop = True
            for gen in self._gen_dfs.keys():
                df = self._gen_dfs[gen]
                if (df is None) or (i >= len(df)):
                    continue
                should_stop = False
                row = df.iloc[i, :].to_numpy()
                if (
                    gen >= 120
                ):  # Note: hack to support latinos genenrations of 121, 122, 123
                    row[0] = row[0][:-2] + str(gen - 120) + row[0][-2:]
                rows.append(row)
            i += 1

        # load next batch
        self._load_dfs()

        return rows

    def get_headers(self):
        pred_header = []
        label_header = []

        for chrom, window_size in sorted(self._window_info.items(), key=lambda x: x[0]):
            pred_header += [
                f"c{chrom}-w{w}-p{p}"
                for w in range(window_size)
                for p in range(self._n_populations)
            ]
            label_header += [f"c{chrom}-w{w}" for w in range(window_size)]

        return ["sample_id"] + pred_header + label_header

    def _load_dfs(self):
        self._finished = True
        for gen in self._gen_df_itrs.keys():
            self._gen_dfs[gen] = None
            sample_ids = []
            preds = []
            labels = []

            window_df_itrs = self._gen_df_itrs[gen]
            for window_df_itr in window_df_itrs:
                try:
                    df = next(window_df_itr)
                except StopIteration:
                    break

                sample_ids.append(df.sample_id.to_numpy())
                preds.append(df.iloc[:, 1:-1].to_numpy())
                labels.append(df.label.to_numpy())

            if not sample_ids:
                logging.info(f"   generation {gen} has no records.")
                continue

            self._finished = False  # not finished as long as a generation has records.
            self._n_populations = preds[0].shape[1]

            assert np.all([np.array_equal(sample_ids[0], sids) for sids in sample_ids])
            sample_ids = sample_ids[0]
            preds = np.hstack(preds)
            labels = np.column_stack(labels)
            logging.info(f"   read generation {gen}: {len(sample_ids)} samples.")

            self._gen_dfs[gen] = pd.concat(
                [
                    pd.DataFrame({"sample_id": sample_ids}),
                    pd.DataFrame(preds),
                    pd.DataFrame(labels),
                ],
                axis=1,
            )

    @staticmethod
    def get_sub_folders(root_dir: Path, prefix):
        dirs = [
            (int(path.name[len(prefix) :]), path)
            for path in root_dir.glob(f"{prefix}*")
        ]
        return sorted(dirs, key=lambda x: x[0])

    def _get_generation_dataframe_iterators(self, input_dir: Path):
        logging.info(f"Input: {input_dir}")
        logging.info(os.listdir(input_dir))

        chrom_dirs = self.get_sub_folders(input_dir, "chr")
        logging.info(f"{chrom_dirs}")

        # generate list of pandas dataframe iterators for each generation.
        # dataframes will be ordered by chromosomes and windows
        gen_df_itrs = {}
        self._window_info = {}
        for chrom, chrom_dir in chrom_dirs:
            logging.info(f"Chromosome {chrom}")
            gen_dirs = self.get_sub_folders(chrom_dir, "gen")
            logging.info(os.listdir(chrom_dir))
            n_generations = len(gen_dirs)
            logging.info(f"Number of generations in {chrom} is {n_generations}")

            for gen, gen_dir in gen_dirs:
                window_dirs = self.get_sub_folders(gen_dir, "window-")
                n_windows = len(window_dirs)
                if chrom not in self._window_info:
                    self._window_info[chrom] = n_windows
                else:
                    assert self._window_info[chrom] == n_windows

                if gen not in gen_df_itrs:
                    gen_df_itrs[gen] = []

                for _, window_dir in window_dirs[:n_windows]:
                    pred_file = window_dir / "prediction.tsv.gz"
                    df_itr = pd.read_csv(pred_file, sep="\t", chunksize=1024)
                    gen_df_itrs[gen].append(df_itr)

        return gen_df_itrs


def merge_base_layer(input_dir: Path, output_dir: Path):
    gen_df_itrs = GenerationDataFrameIterator(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_rows = 0

    with gzip.open(output_dir / "predictions.tsv.gz", "wt") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(gen_df_itrs.get_headers())
        while rows := gen_df_itrs.get_next_batch():
            writer.writerows(rows)
            n_rows += len(rows)
            logging.info(f"{n_rows} written.")

    with open(output_dir / "parameters.json", mode="wt") as fout:
        parameters = {
            "windows": gen_df_itrs.window_info,
            "n-populations": gen_df_itrs.n_populations,
        }
        json.dump(parameters, fout, indent=4)
