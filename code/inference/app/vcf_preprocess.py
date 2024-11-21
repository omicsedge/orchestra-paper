import logging
import typing
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Dict, List, Tuple, cast

import dask.dataframe as dd
import numpy as np
import pandas as pd
from cyvcf2 import VCF, Variant
from dask import delayed

logging.getLogger().setLevel(logging.INFO)


Genotype = Tuple[int, int, bool]  # in fact it's a list but that's just a poor design
# choice from the authors of the lib and it doesn't matter here. Saying it's a tuple makes
# it possible have better type narrowing.


@dataclass(frozen=True)
class Snp:
    """Represent SNP."""

    chrom: int
    pos: int
    ref: str
    alt: str

    def __repr__(self) -> str:
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"

    @classmethod
    def from_variant(cls, v: Variant):
        chrom = v.CHROM[len("chr") :] if "chr" in v.CHROM else v.CHROM
        return cls(int(chrom), v.POS, v.REF, v.ALT[0] if len(v.ALT) else v.ALT)

    @classmethod
    def from_row(cls, row: pd.Series):
        return cls(int(row["chr"]), int(row["pos"]), row["ref"], row["alt"])

    def __gt__(self, other):
        return self.chrom > other.chrom or (
            self.chrom == other.chrom and self.pos > other.pos
        )


class VcfWindowizer:
    """Convert input file into windows."""

    def __init__(self, model_dir: Path) -> None:
        self._model_dir = model_dir
        self._df_snps = self._load_model_snps()

    def get_chromosomes(self) -> typing.List[int]:
        return list(self._df_snps.chr.unique())

    def get_windows(self, chromosome: str) -> typing.List[int]:
        return list(self._df_snps[self._df_snps.chr == chromosome].window.unique())

    def _load_model_snps(self) -> pd.DataFrame:
        snp_files = self._model_dir.rglob("chr*/window-*/snps.tsv")
        dfs = [delayed(self._read_window(snp_file)) for snp_file in snp_files]
        df = dd.from_delayed(dfs).compute()
        return df.sort_values(["chr", "window", "pos"])

    @classmethod
    def _read_window(cls, filename) -> pd.DataFrame:
        # reads each csv file to a pandas.DataFrame
        df = pd.read_csv(filename, sep="\t")
        df["window"] = int(filename.parts[-2][len("window-") :])
        return df

    def read(self, filepath: Path, chromosome: int):
        """
        Read VCF File and convert it to windowed genotype data.

        Args
            filepath: VCF file path (full or relative)
            chromosome: chromosome to read

        Returns
            dictionary with (chromosome, window) as key and numpy array with
            the shape (n_samples, n_window_size) as value.
        """
        vcf = VCF(filepath / f"chr{chromosome}.vcf.gz")
        extracted: Dict[Tuple[int, str], List[List[Tuple[int, int]]]] = defaultdict(
            list
        )
        for _, row in self._df_snps[self._df_snps.chr == chromosome].iterrows():
            snp_model = Snp.from_row(row)
            window: str = row["window"]

            while True:
                variant = next(vcf)
                snp_vcf = Snp.from_variant(variant)

                if snp_vcf > snp_model:
                    raise ValueError(f"SNP {snp_model} not available")

                if snp_vcf == snp_model:
                    gt: List[Tuple[int, int]] = [
                        _gt[:2] for _gt in cast(List[Genotype], variant.genotypes)
                    ]
                    extracted[(chromosome, window)].append(gt)
                    break
        concatenated_array = np.concatenate([*extracted.values()])
        data = concatenated_array.reshape(concatenated_array.shape[0], -1).T

        samples = [sample for sample in product(vcf.samples, [0, 1])]
        return samples, data


def _read_genotypes(
    windowizer: VcfWindowizer, chromosomes: List[int], input_prefix: str
) -> dict:
    data = []
    samples = None
    for chrom in chromosomes:
        samples_, data_ = windowizer.read(input_prefix, chrom)

        # get samples
        if samples is None:
            samples = samples_
        else:
            assert samples == samples_

        # merge data
        data.append(data_)

    return zip(samples, np.concatenate(data, axis=1).tolist())


def vcf_preprocess(
    input_dir: Path,
    model_dir: Path,
    output_dir: Path,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    # init windowizer
    windowizer = VcfWindowizer(model_dir)

    chroms = windowizer.get_chromosomes()

    for sample_batch in input_dir.iterdir():
        logging.info(f"extracting samples in {sample_batch}")

        genotypes = _read_genotypes(windowizer, chroms, sample_batch)

        logging.info(f"writing to {output_dir / f'{sample_batch.parts[-1]}.tsv'}")

        with open(output_dir / f"{sample_batch.parts[-1]}.tsv", "wb") as fout:
            for line, haplotype in genotypes:
                fout.write(
                    (
                        "\t".join(
                            map(
                                str,
                                [
                                    "sample_list_1",
                                    "sample_list_1_file",
                                    *line,
                                    *haplotype,
                                ],
                            )
                        )
                        + "\n"
                    ).encode()
                )
