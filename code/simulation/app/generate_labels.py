import csv
import gzip

import click
import numpy as np
import pandas as pd
import tskit


@click.command()
@click.option(
    "--output-tree-file",
    required=True,
    help="Path to the output tree file.",
)
@click.option(
    "--monomorphic-vts-file",
    required=True,
    help="Path to the monomorphic VTS file.",
)
@click.option(
    "--polymorphic-vts-file",
    required=True,
    help="Path to the polymorphic VTS file.",
)
@click.option(
    "--samples-order-file",
    required=True,
    help="Path to the samples order file.",
)
@click.option(
    "--output-labels-file",
    required=True,
    help="Path to the output labels file.",
)
@click.option(
    "--cohort",
    required=True,
    help="Dataset cohort: test/train usually.",
)
@click.option(
    "--gen",
    required=True,
    help="Simulation generation.",
)
@click.option(
    "--initial-genomes",
    required=True,
    help="Former NUMBER_OF_INITIAL_GENOMES",
    type=int,
)
@click.option(
    "--total-genomes",
    required=True,
    help="Former NUMBER_OF_TOTAL_GENOMES",
    type=int,
)
def generate_labels(
    output_tree_file: str,
    monomorphic_vts_file: str,
    polymorphic_vts_file: str,
    samples_order_file: str,
    output_labels_file: str,
    cohort: str,
    gen: str,
    initial_genomes: int,
    total_genomes: int,
) -> None:
    """Generate labels for genetic variants.

    Args:
        output_tree_file: Path to the tree sequence file
        monomorphic_vts_file: Path to monomorphic variants file
        polymorphic_vts_file: Path to polymorphic variants file
        samples_order_file: Path to samples order file
        output_labels_file: Path to output labels file
        cohort: Dataset cohort (test/train)
        gen: Simulation generation
        initial_genomes: Number of initial genomes
        total_genomes: Total number of genomes
    """
    # Load and process tree sequence
    ts = _load_tree_sequence(output_tree_file)

    # Load and merge variant data
    df = _load_variant_data(monomorphic_vts_file, polymorphic_vts_file)
    snps = df.to_numpy()

    # Generate sample mappings
    idx_to_pop = _load_sample_mappings(samples_order_file)
    idx_to_pop_func = np.vectorize(idx_to_pop.get)

    # Process samples
    thisgen_samples = ts.samples()[initial_genomes:total_genomes]

    _write_labels(
        output_labels_file,
        ts,
        snps,
        thisgen_samples,
        idx_to_pop_func,
        initial_genomes,
        cohort,
        gen,
        df,
    )


def _load_tree_sequence(tree_file: str) -> tskit.TreeSequence:
    """Load and simplify tree sequence."""
    ts = tskit.load(tree_file)
    return ts.simplify(reduce_to_site_topology=True)


def _load_variant_data(mono_file: str, poly_file: str) -> pd.DataFrame:
    """Load and merge monomorphic and polymorphic variant data."""
    monomorphic = pd.read_csv(
        mono_file,
        sep="\t",
        usecols=[0, 1],
        names=["#CHROM", "POS"],
    )
    monomorphic["type"] = "MONO"

    polymorphic = pd.read_csv(poly_file, sep="\t")
    polymorphic["type"] = "POLY"

    return pd.concat([polymorphic, monomorphic]).sort_values(["#CHROM", "POS"])


def _load_sample_mappings(samples_file: str) -> dict:
    """Load sample ID to population mappings."""
    idx_to_pop = {}
    with open(samples_file, "rt") as fin:
        for line in fin:
            idx, _, pop = line.split()
            idx_to_pop[int(idx)] = pop
    return idx_to_pop


def _write_labels(
    output_file: str,
    ts: tskit.TreeSequence,
    snps: np.ndarray,
    thisgen_samples: np.ndarray,
    idx_to_pop_func: np.vectorize,
    initial_genomes: int,
    cohort: str,
    gen: str,
    df: pd.DataFrame,
) -> None:
    """Write labels to output file."""
    with gzip.open(output_file, "wt") as fpop_out:
        pop_writer = csv.writer(fpop_out, delimiter="\t", lineterminator="\n")

        # Write header
        header = ["#CHROM", "POS"] + [
            f"i{i}_{cohort}_generation{gen}" for i in range(len(thisgen_samples) // 2)
        ]
        pop_writer.writerow(header)

        i_snp = 0
        for tree in ts.trees():
            roots = _get_root_ancestors(tree, thisgen_samples, initial_genomes)
            pop_line = _format_population_line(roots, idx_to_pop_func)

            i_snp = _write_snp_data(pop_writer, snps, i_snp, tree.num_sites, pop_line)

        # Write remaining monomorphic SNPs
        while i_snp < len(df):
            assert snps[i_snp, 2] == "MONO"
            pop_writer.writerow(snps[i_snp, :2].tolist() + pop_line.tolist())
            i_snp += 1


def _get_root_ancestors(
    tree: tskit.Tree, samples: np.ndarray, initial_genomes: int
) -> list:
    """Get root ancestors for given samples."""
    roots = []
    for u in samples:
        while tree.parent(u) >= initial_genomes:
            u = tree.parent(u)
        roots.append(tree.parent(u))
    return roots


def _format_population_line(roots: list, idx_to_pop_func: np.vectorize) -> np.ndarray:
    """Format population line for output."""
    pop_line = idx_to_pop_func(np.array(roots))
    return np.apply_along_axis(
        lambda x: f"{x[0]}|{x[1]}",
        0,
        np.vstack((pop_line[::2], pop_line[1::2])),
    )


def _write_snp_data(
    writer: csv.writer,
    snps: np.ndarray,
    i_snp: int,
    num_sites: int,
    pop_line: np.ndarray,
) -> int:
    """Write SNP data and return updated SNP index."""
    for _ in range(num_sites):
        while snps[i_snp, 2] == "MONO":
            writer.writerow(snps[i_snp, :2].tolist() + pop_line.tolist())
            i_snp += 1

        writer.writerow(snps[i_snp, :2].tolist() + pop_line.tolist())
        i_snp += 1

    return i_snp


if __name__ == "__main__":
    generate_labels()
