import logging
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List

import click
from scripts.split_dataset import split_dataset


def process_chromosome(chr: int, input_file: str, temp_path: Path) -> None:
    chr_dir = temp_path / "chrs"
    chr_dir.mkdir(exist_ok=True)
    output_file = chr_dir / f"chr{chr}.vcf.gz"

    # Run bcftools view
    os.system(f"bcftools view -Oz -r {chr} {input_file} -o {output_file}")

    # Run bcftools index
    os.system(f"bcftools index {output_file}")


@click.command()
@click.option("--start-chr", "-sc", type=int, help="Start chromosome")
@click.option("--end-chr", "-ec", type=str, help="End chromosome")
@click.option("--source-panel", type=str, help="Source panel path")
@click.option("--sample-map", type=str, help="Sample map file")
@click.option("--description", type=str, help="Description")
@click.option("--type", type=str, help="Type", default="random")
@click.option("--n-times", type=int, help="Number of times", default=2)
@click.option("--output-dir", type=Path, help="Output directory", default=".")
def simulate(
    start_chr: int,
    end_chr: int,
    source_panel: str,
    sample_map: str,
    description: str,
    type: str,
    n_times: int,
    output_dir: str,
) -> None:
    chromosomes: List[int] = list(range(start_chr, end_chr + 1))

    logging.info(f"Creating index for source panel {source_panel}")
    os.system(f"bcftools index -f {source_panel}")

    # Create a temporary directory that will be automatically cleaned up
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)

        # Process chromosomes in parallel
        with ThreadPoolExecutor() as executor:
            futures = [
                executor.submit(process_chromosome, chr, source_panel, output_dir)
                for chr in chromosomes
            ]
            # Wait for all tasks to complete
            for future in futures:
                future.result()

        # cohorts = split_dataset(sample_map, description, temp_path)
        # for cohort, sample_map_path in cohorts.items():
        #     # Generate pedigree
        #     os.system(f"bash scripts/generate_pedigree.sh {cohort} {sample_map_path} {temp_path / 'chrs' / f'chr{end_chr}.vcf.gz'} {output_dir}")
            
        #     for chr in chromosomes:
        #         # Extract founders
        #         os.system(f"bash scripts/extract_founders.sh {cohort} {temp_path / 'chrs' / f'chr{chr}.vcf.gz'} {output_dir}")

        #         # Simulate generations
        #         for generation in range(1, 7):
        #             os.system(f"bash scripts/simulate_generation.sh {cohort} {chr} {generation} {type} {n_times}")

