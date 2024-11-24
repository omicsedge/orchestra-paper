import csv
import logging
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List, Tuple

import click

from .app.split_dataset import split_dataset

order = [
    "sample_id",
    "level_3_code",
    "superpopulation_code",
    "level_1_name",
    "level_2_name",
    "level_3_name",
    "data_set",
    "meta_1",
]


def process_chromosome(chr: int, input_file: str, temp_path: Path) -> Tuple[int, Path]:
    chr_dir = temp_path / "chrs"
    chr_dir.mkdir(exist_ok=True)
    output_file = chr_dir / f"chr{chr}.vcf.gz"

    # Run bcftools view
    run_command(f"bcftools view -Oz -r {chr} {input_file} -o {output_file}")

    # Run bcftools index
    run_command(f"bcftools index {output_file}")

    return chr, output_file


def run_command(cmd: str) -> None:
    exit_code = os.system(cmd)
    if exit_code != 0:
        raise RuntimeError(f"Command failed with exit code {exit_code}: {cmd}")


@click.command()
@click.option("--start-chr", "-sc", type=int, help="Start chromosome")
@click.option("--end-chr", "-ec", type=int, help="End chromosome")
@click.option("--source-panel", "-sp", type=str, help="Source panel path")
@click.option("--sample-map", "-sm", type=str, help="Sample map file")
@click.option("--version", "-v", type=str, help="Version")
@click.option("--type", "-t", type=str, help="Type", default="random")
@click.option("--n-times", "-nt", type=int, help="Number of times", default=2)
@click.option(
    "--output-dir",
    "-o",
    type=Path,
    help="Output directory",
    default="/results/simulation",
)
def simulation(
    start_chr: int,
    end_chr: int,
    source_panel: str,
    sample_map: str,
    version: str,
    type: str,
    n_times: int,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    chromosomes: List[int] = list(range(start_chr, end_chr + 1))

    # Create a temporary directory that will be automatically cleaned up
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        chr_files = {}

        with open(sample_map) as f1, open(temp_path / "sample_map.tsv", "w") as f2:
            reader = csv.DictReader(f1, delimiter="\t")
            writer = csv.writer(f2, delimiter="\t")
            writer.writerow(order)

            for row in reader:
                new_row = {
                    "sample_id": row["SampleID"],
                    "level_3_code": row["level1"],
                    "superpopulation_code": row["level1"],
                    "level_1_name": row["level3"],
                    "level_2_name": row["level3"],
                    "level_3_name": row["level3"],
                    "data_set": row["Dataset"],
                    "meta_1": row["Population"],
                }
                data = []
                for i in order:
                    data.append(new_row[i])
                writer.writerow(data)

        sample_map = temp_path / "sample_map.tsv"

        # Process chromosomes in parallel
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_chromosome, chr, source_panel, output_dir)
                for chr in chromosomes
            ]
            # Wait for all tasks to complete
            for future in futures:
                chr, output_file = future.result()
                chr_files[chr] = output_file

        cohorts = split_dataset(sample_map, version, temp_path, output_dir)
        for cohort, sample_map_path in cohorts.items():
            # Generate pedigree
            run_command(
                "bash /code/simulation/scripts/generate_pedigree.sh "
                f"{cohort} {sample_map_path} {chr_files[end_chr]} "
                f"{output_dir} {end_chr} {type} {n_times}"
            )

            for chr in chromosomes:
                # Extract founders
                run_command(
                    "bash /code/simulation/scripts/extract_founders.sh "
                    f"{cohort} {sample_map_path} {chr_files[chr]} {output_dir} {chr}"
                )

                # Simulate generations in parallel
                with ThreadPoolExecutor() as executor:
                    futures = [
                        executor.submit(
                            run_command,
                            f"bash code/simulation/scripts/simulate_generation.sh {cohort} {sample_map_path} {chr_files[chr]} {output_dir} {chr} {generation} {type} {n_times}",
                        )
                        for generation in range(1, 7)
                    ]
                    # Wait for all generations to complete
                    for future in futures:
                        future.result()


if __name__ == "__main__":
    simulation()
