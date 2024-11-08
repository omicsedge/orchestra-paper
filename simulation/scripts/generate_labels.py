import csv
import gzip

import click
import numpy as np
import pandas as pd
import tskit


@click.command()
@click.option(
    "--output-tree-file",
    help="Path to the output tree file.",
)
@click.option(
    "--monomorphic-vts-file",
    help="Path to the monomorphic VTS file.",
)
@click.option(
    "--polymorphic-vts-file",
    help="Path to the polymorphic VTS file.",
)
@click.option(
    "--samples-order-file",
    help="Path to the samples order file.",
)
@click.option(
    "--output-labels-file",
    help="Path to the output labels file.",
)
@click.option(
    "--cohort",
    help="Dataset cohort: test/train usually.",
)
@click.option(
    "--gen",
    help="Simulation generation.",
)
@click.option(
    "--initial-genomes",
    help="Former NUMBER_OF_INITIAL_GENOMES",
    type=int,
)
@click.option(
    "--total-genomes",
    help="Former NUMBER_OF_TOTAL_GENOMES",
    type=int,
)
def run(
    output_tree_file: str,
    monomorphic_vts_file: str,
    polymorphic_vts_file: str,
    samples_order_file: str,
    output_labels_file: str,
    cohort: str,
    gen: str,
    initial_genomes: int,
    total_genomes: int,
):
    # Load data and remove regions without variants
    ts = tskit.load(output_tree_file)
    ts = ts.simplify(reduce_to_site_topology=True)

    # Read monomorphic file
    monomorphic = pd.read_csv(
        monomorphic_vts_file,
        sep="\t",
        usecols=[0, 1],
        names=["#CHROM", "POS"],
    )
    monomorphic["type"] = "MONO"

    # Read polymorfic file
    polymorphic = pd.read_csv(polymorphic_vts_file, sep="\t")
    polymorphic["type"] = "POLY"

    # Merge files
    df = pd.concat([polymorphic, monomorphic]).sort_values(["#CHROM", "POS"])
    snps = df.to_numpy()

    # Generate a mapping between tskit index and sample_id and population
    idx_to_sid = {}
    idx_to_pop = {}
    with open(samples_order_file, "rt") as fin:
        for line in fin:
            idx, sample_id, pop = line.split()
            idx = int(idx)
            idx_to_sid[idx] = sample_id + "." + str(idx % 2)
            idx_to_pop[idx] = pop

    # convert mapping to numpy vectorizer
    # idx_to_sid_func = np.vectorize(idx_to_sid.get)
    idx_to_pop_func = np.vectorize(idx_to_pop.get)

    # Approach:
    # Move along the tree sequence (iterates over tree objects)
    # We check the root of the tree starting at the samples (final generation),
    # so we can track local ancestry t1_start = process_time()
    thisgen_samples = ts.samples()[initial_genomes:total_genomes]

    with gzip.open(output_labels_file, "wt") as fpop_out:
        pop_writer = csv.writer(fpop_out, delimiter="\t", lineterminator="\n")

        header = ["#CHROM", "POS"] + [
            f"i{i}_{cohort}_generation{gen}" for i in range(len(thisgen_samples) // 2)
        ]
        pop_writer.writerow(header)

        i_snp = 0
        for tree in ts.trees():
            # find super parents at root of the tree
            roots = []
            for u in thisgen_samples:
                while tree.parent(u) >= initial_genomes:
                    u = tree.parent(u)
                roots.append(tree.parent(u))

            # write to files
            pop_line = idx_to_pop_func(np.array(roots))
            pop_line = np.apply_along_axis(
                lambda x: f"{x[0]}|{x[1]}",
                0,
                np.vstack((pop_line[::2], pop_line[1::2])),
            )
            # sid_line_ = idx_to_sid_func(roots)
            # sid_line = []
            # for a, b in zip(sid_line_[::2], sid_line_[1::2]):
            #     sid_line.append(a + '|' + b)

            for i in range(tree.num_sites):
                # Write the same line to Monomorphic SNPs.
                while snps[i_snp, 2] == "MONO":
                    pop_writer.writerow(snps[i_snp, :2].tolist() + pop_line.tolist())
                    i_snp += 1

                # Write polymorphic SNPs
                pop_writer.writerow(snps[i_snp, :2].tolist() + pop_line.tolist())
                i_snp += 1

        # Write monomorphic SNPs if left.
        while i_snp < len(df):
            assert snps[i_snp, 2] == "MONO"
            pop_writer.writerow(snps[i_snp, :2].tolist() + pop_line.tolist())
            i_snp += 1


if __name__ == "__main__":
    run()
