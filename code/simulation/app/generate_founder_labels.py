import csv
import gzip
import subprocess
from typing import Dict, List, Tuple

import click


def get_samples(vcf_file: str) -> List[str]:
    """Get sample IDs from VCF file.

    Args:
        vcf_file: Path to the VCF file

    Returns:
        List of sample IDs found in the VCF file
    """
    results = subprocess.run(
        ["bcftools", "query", "-l", vcf_file], stdout=subprocess.PIPE, check=True
    )
    return results.stdout.decode().strip().split("\n")


def get_loci(vcf_file: str) -> List[List[str]]:
    """Get chromosome and position information from VCF file.

    Args:
        vcf_file: Path to the VCF file

    Returns:
        List of [chromosome, position] pairs
    """
    results = subprocess.run(
        ["bcftools", "query", "-f", "%CHROM\t%POS\n", vcf_file],
        stdout=subprocess.PIPE,
        check=True,
    )
    return [loci.split("\t") for loci in results.stdout.decode().strip().split("\n")]


def get_sample_map(sample_map_file: str) -> Dict[str, str]:
    """Load sample mapping from file into a dictionary.

    Args:
        sample_map_file: Path to the sample mapping file

    Returns:
        Dictionary mapping sample IDs to their level 3 codes
    """
    with open(sample_map_file) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        sample_id_col = header.index("sample_id")
        level_3_code_col = header.index("level_3_code")

        return {line[sample_id_col]: line[level_3_code_col] for line in reader}


@click.command()
@click.option("--vcf-file", required=True, help="Path to the VCF file")
@click.option(
    "--sample-map-file", required=True, help="Path to the sample mapping file"
)
@click.option("--output-file", required=True, help="Path to the output gzip file")
def generate_founder_labels(
    vcf_file: str, sample_map_file: str, output_file: str
) -> None:
    """Generate founder labels from VCF and sample mapping data.

    Processes VCF file and sample mapping to create a gzipped output file containing
    founder labels for each chromosome position.
    """
    sample_ids = get_samples(vcf_file)
    loci = get_loci(vcf_file)
    sample_map = get_sample_map(sample_map_file)

    click.echo("Generating founder labels...")

    with gzip.open(output_file, mode="wt") as fout:
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

        # Filter for founders (samples that exist in the mapping)
        founders = [_id for _id in sample_ids if _id in sample_map]

        # Create header row
        writer.writerow(["#CHROM", "POS"] + founders)

        # Create founder labels (doubled because of diploid organisms)
        founder_labels = [
            "|".join([sample_map[sample_id]] * 2) for sample_id in founders
        ]

        # Write data rows
        for chrom, pos in loci:
            writer.writerow([chrom, pos] + founder_labels)

    click.echo(f"Successfully wrote founder labels to {output_file}")


if __name__ == "__main__":
    generate_founder_labels()
