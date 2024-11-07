import csv
import gzip
import subprocess

import click


def get_samples(vcf_file):
    """Get sample ids from VCF"""
    results = subprocess.run(
        ["bcftools", "query", "-l", vcf_file], stdout=subprocess.PIPE
    )
    samples = results.stdout.decode().strip().split("\n")
    return samples


def get_loci(vcf_file):
    """Get chrom and pos."""
    results = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\n", vcf_file], stdout=subprocess.PIPE
    )
    return [loci.split("\t") for loci in results.stdout.decode().strip().split("\n")]


def get_sample_map(sample_map_file):
    """Get sample map as an dictionary."""
    with open(sample_map_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        sample_id_col = header.index("sample_id")
        level_3_code_col = header.index("level_3_code")

        sample_map = {line[sample_id_col]: line[level_3_code_col] for line in reader}

    return sample_map


@click.command()
@click.option("--vcf-file", required=True, help="VCF File.")
@click.option("--sample-map-file", required=True, help="Sample Map File.")
@click.option("--output-file", required=True, help="Output gzip file")
def get_founders(vcf_file, sample_map_file, output_file):
    """Fix labels in the label file."""
    sample_ids = get_samples(vcf_file)
    loci = get_loci(vcf_file)
    sample_map = get_sample_map(sample_map_file)

    click.echo("Fixing label file.")

    with gzip.open(output_file, mode="wt") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        # find missing sample_ids. only founders are missing so take the
        founders = [_id for _id in sample_ids if _id in sample_map]

        # write new header
        new_header = ["#CHROM", "POS"] + founders
        writer.writerow(new_header)

        # missing sample label list
        founder_labels = [
            "|".join([sample_map[sample_id]] * 2) for sample_id in founders
        ]

        for chrom, pos in loci:
            writer.writerow([chrom, pos] + founder_labels)


if __name__ == "__main__":
    get_founders()
