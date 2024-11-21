import click


@click.command()
@click.option(
    "--vcf-file",
    required=True,
    help="VCF File.",
    default="tools/simulation/results/chr20/test/test_set.vcf",
)
@click.option(
    "--remove-variants-file",
    required=True,
    help="File with variants to remove.",
    default="tools/simulation/results/chr20/test/remove_variants.txt",
)
@click.option(
    "--output-file",
    required=True,
    help="Output file",
    default="tools/simulation/output.vcf",
)
def run(vcf_file, remove_variants_file, output_file):
    """Fix labels in the label file."""
    click.echo("Read remove variants.")
    with open(remove_variants_file, mode="rt") as fin:
        remvars = [line.strip() for line in fin.readlines()]

    click.echo("Process file.")
    with open(vcf_file, mode="rt") as fin, open(output_file, mode="wt") as fout:
        for i, line in enumerate(fin.readlines()):
            # write comments
            if line.startswith("#"):
                fout.write(line)
                continue

            if i % 100000 == 0:
                click.echo(f"  {i} lines")

            # line iterator
            iter_line = iter(line)

            # find first tab (i.e. skip chrom)
            while next(iter_line) != "\t":
                pass

            # read pos
            pos = ""
            while (ch := next(iter_line)) != "\t":
                pos += ch

            if not pos in remvars:
                fout.write(line)


if __name__ == "__main__":
    run()
