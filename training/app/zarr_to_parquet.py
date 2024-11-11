""""""

import logging
import os
from collections import defaultdict
from functools import partial
from itertools import product
from pathlib import Path
from typing import Callable, Dict, List

import numpy as np
import pandas as pd
import xarray
from joblib import Parallel, delayed

logging.basicConfig(level=logging.INFO)


class PreProcess:
    """
    Handles preprocessing of genetic data by:
    1. Loading and decompressing ZARR files
    2. Converting data into windowed segments
    3. Saving processed data as parquet files
    """

    def __init__(
        self,
        input_dir: Path,
        output_dir: Path,
        chromosome: int,
        window_size: int,
        mapper_fn: Callable[[np.ndarray], np.ndarray],
    ) -> None:
        """
        Args:
            input_dir: Directory containing input ZARR files
            output_dir: Directory to save processed data
            chromosome: Chromosome number to process
            window_size: Size of windows to segment data
            mapper_fn: Function to map between population levels
        """
        self._input_dir = input_dir
        self._output_dir = output_dir
        self._chromosome = chromosome
        self._window_size = window_size
        self._mapper_fn = mapper_fn

    @staticmethod
    def get_level_mapper(
        pop_map: str, level: int
    ) -> Callable[[np.ndarray], np.ndarray]:
        """
        Creates a mapping function to convert between population hierarchy levels.
        Args:
            pop_map: Path to TSV file containing population hierarchy mappings
            level: The level label.
        Returns:
            Vectorized mapping function to convert population IDs
        """
        col_name: str = f"level_{level}_name"
        df = pd.read_csv(pop_map, sep="\t")
        level_to_id = dict(
            {p: i for i, p in enumerate(df[col_name].drop_duplicates().to_list())}
        )
        level_3_id_to_new_id = df[col_name].map(level_to_id).to_dict()
        return np.vectorize(level_3_id_to_new_id.get)

    def run(self) -> None:
        zarr_files = self._get_zarr_files()
        self._decompress_zarr_files(zarr_files)

        # convert each generation to windowed CSVs
        for g in zarr_files.keys():
            files = [f[1] for f in sorted(zarr_files[g].items(), key=lambda x: x[0])]
            self._preprocess_generation(files, g)

    def _get_zarr_files(self) -> Dict[int, Dict[int, str]]:
        """
        Scans input directory for ZARR files and organizes them by generation.
        Returns:
            Dictionary mapping generation numbers to ZARR file paths
        """
        zarr_files = defaultdict(dict)
        gens = self._input_dir.glob("gen*")

        for gen in gens:
            g = int(gen.name[len("gen") :])

            zarr_file = list(gen.glob(f"chr{self._chromosome}.zarr.tar.gz"))[0]
            zarr_files[g][0] = str(zarr_file)[: -len(".tar.gz")]

        logging.info(zarr_files)

        return zarr_files

    def _preprocess_generation(self, zarr_files: List[str], generation: int) -> None:
        """
        Processes all data for a single generation:
        1. Loads genetic data from ZARR files
        2. Extracts sample information, genotypes, and variant data
        3. Segments data into windows
        4. Processes each window in parallel

        Args:
            zarr_files: List of ZARR files for this generation
            generation: Generation number being processed
        """
        logging.info(f"Pre-processing generation: {generation}")
        dss = []
        for zfile in zarr_files:
            dss.append(xarray.open_zarr(zfile))

        # extract data from zarr file
        sample_names = np.hstack([ds.samples.values for ds in dss])
        genotypes = np.vstack([ds.call_genotype.values.swapaxes(0, 1) for ds in dss])
        labels = np.vstack([ds.admixture.values.swapaxes(0, 1) for ds in dss])
        variant_alleles = np.vstack(
            [
                (
                    ds.variant_allele.values[np.newaxis].T
                    if len(ds.variant_allele.values.shape) == 1
                    else ds.variant_allele.values
                )
                for ds in dss
            ]
        )
        variant_position = np.vstack(
            [ds.variant_position.values[np.newaxis].T for ds in dss]
        )

        # logging info
        logging.info(f"Number of samples: {sample_names.shape[0]}")
        logging.info(f"Number of SNPs: {genotypes.shape[1]}")

        # snps: create array of chr, pos, ref and alt
        chromosomes = np.repeat([[self._chromosome]], variant_position.shape[0], axis=0)
        snps = np.concatenate((chromosomes, variant_position, variant_alleles), axis=1)

        # sample_names
        sample_names = pd.Series(
            [
                "_".join(phase_name[::-1])
                for phase_name in product(("0", "1"), sample_names)
            ],
            name="Sample",
        )

        genotypes = np.vstack([genotypes[:, :, 0], genotypes[:, :, 1]])
        labels = np.array(np.vstack([labels[:, :, 0], labels[:, :, 1]]), dtype=np.int32)

        length = genotypes.shape[-1]

        breakpoints = list(range(0, length, self._window_size))
        windows = list(zip(breakpoints, breakpoints[1:] + [length]))

        # remove last window if partial
        if windows[-1][1] - windows[-1][0] != self._window_size:
            windows = windows[:-1]
        logging.info(f"Number of windows (chr{self._chromosome}) :{len(windows)}")

        process_window_func = partial(
            self._preprocess_window,
            generation,
            sample_names,
            genotypes,
            labels,
            snps,
        )
        Parallel(n_jobs=(os.cpu_count() // 4) * 3, prefer="threads")(
            delayed(process_window_func)(i, ws, we)
            for i, (ws, we) in enumerate(windows)
        )

        logging.info(f"Pre-processing completed: {generation}")

    def _preprocess_window(
        self,
        generation: int,
        sample_names: pd.Series,
        genotypes: np.ndarray,
        labels: np.ndarray,
        snps: np.ndarray,
        w: int,
        ws: int,
        we: int,
    ) -> None:
        """
        Processes a single window of genetic data:
        1. Extracts window-specific data
        2. Converts labels to hard assignments
        3. Saves data as parquet file
        4. For generation 0, also saves SNP information

        Args:
            generation: Generation number
            sample_names: Sample identifiers
            genotypes: Genetic data matrix
            labels: Population labels
            snps: SNP information
            w: Window index
            ws: Window start position
            we: Window end position
        """
        logging.info(f"Window {w} (gen{generation}): {ws}-{we}")
        genotypes_ = genotypes[:, ws:we]
        labels_ = labels[:, ws:we]
        snps_ = snps[ws:we]

        labels_hard = np.apply_along_axis(lambda x: np.bincount(x).argmax(), 1, labels_)
        labels_hard = self._mapper_fn(labels_hard)

        data_path = (
            self._output_dir / f"chr{self._chromosome}/gen{generation}/window-{w}"
        )
        data_path.mkdir(parents=True, exist_ok=True)

        df = pd.concat(
            [
                pd.DataFrame(sample_names.to_numpy(), columns=["sample_id"], dtype=str),
                pd.DataFrame({"chromosome": [self._chromosome] * len(sample_names)}),
                pd.DataFrame({"window": [w] * len(sample_names)}),
                pd.DataFrame(
                    genotypes_,
                    dtype=int,
                    columns=[str(i) for i in range(genotypes_.shape[1])],
                ),
                pd.DataFrame(labels_hard, columns=["label"], dtype=int),
            ],
            axis=1,
        )
        df.to_parquet(data_path / "data.parquet")

        if generation == 0:
            pd.DataFrame(snps_).to_csv(
                data_path / "snps.tsv",
                sep="\t",
                header=["chr", "pos", "ref", "alt"],
                index=False,
            )

        logging.info(f"Window {w} (gen{generation}) completed.")

    def _decompress_zarr_files(self, filenames: Dict[int, Dict[int, str]]) -> None:
        """
        Decompresses ZARR files from their tar.gz format.
        Args:
            filenames: Dictionary of compressed ZARR file paths
        """
        for split in filenames.values():
            for filename in split.values():
                path = Path(filename)
                os.system(f"tar -C {path.parent} -xf {filename}.tar.gz")


def zarr_to_parquet(
    input_dir: Path,
    output_dir: Path,
    chromosome: int,
    window_size: int,
    level: int,
    pop_map: str,
):
    mapper_fn = PreProcess.get_level_mapper(pop_map, level)

    PreProcess(
        input_dir=input_dir,
        output_dir=output_dir,
        chromosome=chromosome,
        window_size=window_size,
        mapper_fn=mapper_fn,
    ).run()
