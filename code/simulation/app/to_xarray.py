import csv
import gzip
import os

import click
import dask.dataframe as dd
import numpy as np
import sgkit
import xarray as xr
from sgkit.io.vcf import vcf_to_zarr


def index_data_file(data_file):
    os.system(f"bcftools index -tf {data_file}")


def add_indices(ds):
    # variant_contig or chromosome is not set because this is single choromosome
    ds = ds.set_index({"variants": ("variant_position")})
    ds = ds.set_index({"samples": ("sample_id")})

    # fix data types (remvove since not used)
    try:
        ds = ds.drop_vars("variant_id")
    except ValueError:
        pass
    try:
        ds = ds.drop_vars("sample_id")
    except ValueError:
        pass
    try:
        ds = ds.drop_vars("variant_allele")
    except ValueError:
        pass

    return ds


def clean_ds(ds):
    try:
        ds = ds.drop_vars("variant_id")
    except ValueError:
        pass
    try:
        ds = ds.drop_vars("sample_id")
    except ValueError:
        pass
    try:
        ds = ds.drop_vars("variant_allele")
    except ValueError:
        pass

    # ds["samples"] = ds.samples.astype(str)
    # ds["variant_id"] = ds.variant_id.astype(str)
    # ds["variant_allele"] = ds.variant_allele.astype(str)

    return ds


def add_labels_memef(filename, label_file, population_map):
    # read population map
    df_map = dd.read_csv(population_map, sep="\t").compute()
    code_to_id = (
        df_map[["id", "level_3_code"]].set_index("level_3_code")["id"].to_dict()
    )

    # convert population labels to code
    with gzip.open(label_file, "rt") as fin:
        reader = csv.reader(fin, delimiter="\t")

        # open dataset
        ds = xr.open_zarr(filename)

        # verify samples match
        header = next(reader)
        assert np.all(ds.samples.values == header[2:])

        # placeholder admixture variable
        ds = ds.assign(
            admixture=(
                ["variants", "samples", "ploidy"],
                np.zeros(ds.call_genotype.shape, dtype=np.int8),
            )
        )

        # calculate buffer size
        buffer_lines = int(
            ((16 * 1024 * 1024 * 1024) / 32) // (ds.admixture.shape[1] * 8)
        )
        print(f"Buffer lines: {buffer_lines}")

        # save & close
        ds = clean_ds(ds)
        ds = ds.to_zarr(filename, mode="a")
        ds.close()

        # add admixture labels
        n_line = 0
        buffer = []
        for line in reader:
            tmp = []
            for geno_pop in line[2:]:
                l_pop, r_pop = geno_pop.split("|")
                tmp.append([code_to_id[l_pop], code_to_id[r_pop]])

            buffer.append(tmp)

            n_line += 1
            if n_line % buffer_lines == 0:
                # save buffer
                ds = xr.open_zarr(filename)
                ds.admixture[(n_line - buffer_lines) : n_line, :, :] = np.array(buffer)
                ds = clean_ds(ds)
                ds = ds.to_zarr(filename, mode="a")
                ds.close()
                buffer = []

            if n_line % 1000 == 0:
                print(f"..{n_line}")

        # save remaining buffer
        ds = xr.open_zarr(filename)
        ds.admixture[(n_line - len(buffer)) : n_line, :, :] = np.array(buffer)
        ds = clean_ds(ds)
        ds = ds.to_zarr(filename, mode="a")
        ds.close()


def _sample_id_to_gen(sample_id: str) -> int:
    """Convert sample_id to generation"""
    gen_str = "generation"
    if gen_str not in sample_id:
        return 0  # generation 0 (i.e. founder)

    idx = sample_id.index(gen_str)
    gen = sample_id[idx + len(gen_str) :]

    return int(gen)


def add_generations(ds):
    """Add generation for each sample."""
    generations = np.vectorize(_sample_id_to_gen)(ds.samples.values)
    ds = ds.assign(generation=(["samples"], generations))

    return ds


@click.command()
@click.option("--data-file", help="VCF file with genotype data", required=True)
@click.option(
    "--label-file", help="Label file with population for each SNP.", required=True
)
@click.option(
    "--population-map", help="Population map to convert code to id", required=True
)
@click.option("--output-path", help="Path to store output zarr file.", required=True)
def to_xarray(
    data_file,
    label_file,
    population_map,
    output_path,
):
    """Convert Datafile (VCF) and Label file to XArray format."""
    # add index file for the data VCF
    print("Index data file...")
    index_data_file(data_file)

    # convert data file to zarr
    print("VCF to Zarr conversion...")
    vcf_to_zarr(data_file, output_path, max_alt_alleles=1)

    # load dataset
    print("Load dataset...")
    ds = sgkit.load_dataset(output_path)

    # add indices
    print("Add indices...")
    ds = add_indices(ds)
    ds.to_zarr(output_path, mode="a")

    # add labels
    print("Add labels...")
    add_labels_memef(output_path, label_file, population_map)
    ds = sgkit.load_dataset(output_path)  # reload since adding admixture closes

    # add generations
    print("Add generations...")
    ds = add_generations(ds)
    ds = clean_ds(ds)
    ds.to_zarr(output_path, mode="a")


if __name__ == "__main__":
    to_xarray()
