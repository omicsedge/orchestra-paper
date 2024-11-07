import logging
import os
import tarfile
import typing
from collections import defaultdict
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, cast

import boto3
import dask.dataframe as dd
import numpy as np
import pandas as pd
from cyvcf2 import VCF, Variant
from dask import delayed

logger = logging.getLogger()
logger.setLevel(logging.INFO)


@dataclass
class SampleInfo:
    "Information about Input Sample."
    user_id: str
    file_id: str
    input_prefix: str


# This is the data returned by pyarrow when doing `to_pylist`
DictVariant = typing.TypedDict(
    "DictVariant",
    {
        "chrom": str,
        "pos": int,
        "ref": str,
        "alt": str,
        "gt1": typing.Literal[0, 1],
        "gt2": typing.Literal[0, 1],
    },
)
Genotype = Tuple[int, int, bool]  # in fact it's a list but that's just a poor design
# choice from the authors of the lib and it doesn't matter here. Saying it's a tuple makes
# it possible have better type narrowing.

VARIANTS_CHUNK_SIZE = 300


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


class Windowizer:
    """Convert input file into windows."""

    def __init__(self, model_dir: str) -> None:
        self._model_dir = Path(model_dir)
        self._df_snps = self._load_model_snps()

    def get_chromosomes(self) -> typing.List[int]:
        return list(self._df_snps.chr.unique())

    def get_windows(self, chromosome: str) -> typing.List[int]:
        return list(self._df_snps[self._df_snps.chr == chromosome].window.unique())

    def read(self, filepath: str, chromosome: int) -> pd.DataFrame:
        raise NotImplementedError()

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


class VcfWindowizer(Windowizer):
    def _download_chromosome(self, bucket: str, key: str):
        """Download chromosome vcf to temp file"""
        s3 = boto3.client("s3")
        temp_path = "/tmp/temp.vcf.gz"
        print(f"s3://{bucket}/{key}")
        s3.download_file(bucket, key, temp_path)
        return temp_path

    def read(self, bucket: str, filepath: str, chromosome: int):
        """
        Read VCF File and convert it to windowed genotype data.

        Args
            filepath: VCF file path (full or relative)
            chromosome: chromosome to read

        Returns
            dictionary with (chromosome, window) as key and numpy array with
            the shape (n_samples, n_window_size) as value.
        """
        vcf = VCF(
            self._download_chromosome(
                bucket, os.path.join(filepath, f"chr{chromosome}.vcf.gz")
            )
        )
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


def read_genotypes(
    windowizer: Windowizer, chromosomes: List[int], bucket: str, input_prefix: str
) -> dict:
    data = []
    samples = None
    for chrom in chromosomes:
        samples_, data_ = windowizer.read(
            bucket,
            input_prefix,
            chrom,
        )

        # get samples
        if samples is None:
            samples = samples_
        else:
            assert samples == samples_

        # merge data
        data.append(data_)

    return zip(samples, np.concatenate(data, axis=1).tolist())


def get_sample_info(sample: dict):
    return SampleInfo(
        sample["user_id"],
        sample["genome_file_id"],
        os.path.join(sample["user_id"], "genome-files", sample["genome_file_id"]),
    )


def preprocess(
    input_dir: str,
    model_dir: str,
    output_dir: str,
):
    # extract model
    with open(model_dir + "model.tar.gz", "rb") as model_obj:
        with tarfile.open(fileobj=model_obj, mode="r:gz") as tar_obj:
            tar_obj.extractall(path=model_dir)

    # init windowizer
    windowizer = VcfWindowizer(model_dir)
    chroms = windowizer.get_chromosomes()
    print(f"chromosomes: {chroms}")

    print(os.listdir(model_dir))
    path = Path(input_dir)
    for file in path.glob("*.csv"):
        print(f"extacting samples in {file}")
        for file in path.glob("*.csv"):
            print("extract samples...")
            samples = []
            with open(output_dir + file.name, "wb") as fout, open(file, "rt") as fin:
                for line in fin:
                    data = line.strip().split(",")
                    sample_info = SampleInfo(
                        user_id=data[0],
                        file_id=data[0] + "_file",
                        input_prefix=data[2],
                    )

                    genotypes = read_genotypes(
                        windowizer, chroms, data[1], sample_info.input_prefix
                    )

                    for line, haplotype in genotypes:
                        fout.write(
                            (
                                "\t".join(
                                    map(
                                        str,
                                        [
                                            sample_info.user_id,
                                            sample_info.file_id,
                                            *line,
                                            *haplotype,
                                        ],
                                    )
                                )
                                + "\n"
                            ).encode()
                        )
