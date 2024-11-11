# Packages
import click
import tskit


@click.command()
@click.option(
    "--input-file",
    help="Path to the input file.",
    type=str,
)
@click.option(
    "--output-file",
    help="Path to the output file.",
    type=str,
)
@click.option(
    "--extract-individuals",
    help="Former NUMBER_EXTRACT_INDIVIDUALS",
    type=int,
)
def run(input_file: str, output_file: str, extract_individuals: int):
    # Load data and remove regions without variants
    ts = tskit.load(input_file)
    ts = ts.simplify(reduce_to_site_topology=True)

    # Extract VCF
    with open(output_file, "w") as vcf_file:
        ts.write_vcf(vcf_file, individuals=list(range(0, extract_individuals, 1)))


if __name__ == "__main__":
    run()
