"""
Script for model training
"""

import logging
import os
from pathlib import Path
from typing import List, Tuple

import click
from joblib import Parallel, delayed

from .app import (
    AncestryModel,
    evaluate,
    merge_base_layer,
    predict_window,
    smooth_inference,
    train_base_layer,
    train_smooth_layer,
    zarr_to_parquet,
)


@click.command()
@click.option(
    "--simulated-data-path",
    "-sd",
    help="Simulated dataset path",
    type=Path,
)
@click.option("--start-chr", "-sc", type=int, help="Start chromosome")
@click.option("--end-chr", "-ec", type=int, help="End chromosome")
@click.option(
    "--output-dir",
    "-o",
    help="Output for data",
    type=Path,
    default="/results/training",
)
@click.option(
    "--window-size",
    "-ws",
    help="Window size",
    type=int,
    default=600,
)
@click.option(
    "--level",
    "-l",
    help="Level",
    type=int,
    default=3,
)
@click.option(
    "--version",
    "-v",
    help="Version",
    type=str,
    default="0.01",
)
@click.option(
    "--cohorts",
    "-c",
    multiple=True,
    default=["train", "test"],
    show_default=True,
    help="List of cohorts to process.",
)
@click.option(
    "--epochs",
    "-e",
    help="Epochs",
    type=int,
)
def training(
    simulated_data_path: Path,
    start_chr: int,
    end_chr: int,
    output_dir: Path,
    window_size: int,
    level: int,
    version: str,
    cohorts: List[str],
    epochs: int,
) -> None:
    """
    Train the model using the provided parameters.

    Args:
        simulated_data_path: Path to simulated dataset
        start_chr: Starting chromosome number
        end_chr: Ending chromosome number
        output_dir: Output directory for data
        window_size: Size of the window
        level: Level parameter
        version: Version string
        cohorts: List of cohorts to process
    """
    logging.info("Run preprocessing, converting zarr to parquet + windows")
    logging.info(f"Chromosomes: {start_chr}-{end_chr}")

    # Setup directories and parameters
    output_dir = output_dir / f"m{start_chr}.{end_chr}"
    output_dir.mkdir(parents=True, exist_ok=True)

    chromosomes: List[int] = list(range(start_chr, end_chr + 1))
    cpu_count: int = os.cpu_count()
    dataset: Path = output_dir / f"w{window_size}_parquet"

    # Step 1: Convert Zarr to Parquet
    _convert_zarr_to_parquet(
        simulated_data_path,
        chromosomes,
        cohorts,
        dataset,
        window_size,
        level,
        cpu_count,
    )

    # Step 2: Train base model
    train_base_layer(
        input_dir=dataset / "train", output_dir=output_dir / "base_layer_model"
    )

    # Step 3: Predict base layer
    _predict_base_layer(dataset, output_dir, chromosomes, cohorts, cpu_count)

    # Step 4: Merge base predictions
    for cohort in cohorts:
        merge_base_layer(
            input_dir=output_dir / "base_inference" / cohort,
            output_dir=output_dir / "base_inference" / cohort,
        )

    # Step 5: Train smoothing model
    train_smooth_layer(
        pop_map=simulated_data_path / "population_map.tsv",
        training_dir=output_dir / "base_inference" / "train",
        validation_dir=output_dir / "base_inference" / "test",
        output_dir=output_dir / "smooth_layer_model",
        epochs=epochs,
        level=level,
    )

    # Step 6: Predict smooth layer for test data
    smooth_inference(
        smooth_layer_model=output_dir / "smooth_layer_model",
        base_layer_pred=output_dir / "base_inference" / "test" / "predictions.tsv.gz",
        base_layer_params=output_dir / "base_inference" / "test" / "parameters.json",
        output_dir=output_dir / "smooth_inference_test",
    )

    # Step 7: Save model
    _save_model(simulated_data_path, output_dir, version, chromosomes)

    # Step 8: Evaluate results
    # _evaluate_results(output_dir)


def _convert_zarr_to_parquet(
    simulated_data_path: Path,
    chromosomes: List[int],
    cohorts: List[str],
    dataset: Path,
    window_size: int,
    level: int,
    cpu_count: int,
) -> None:
    """Helper function to convert Zarr files to Parquet format."""
    Parallel(n_jobs=cpu_count, prefer="threads")(
        delayed(zarr_to_parquet)(
            pop_map=str(simulated_data_path / "population_map.tsv"),
            input_dir=simulated_data_path / cohort / "zarr-files" / f"chr{chrom}",
            output_dir=dataset / cohort,
            chromosome=chrom,
            window_size=window_size,
            level=level,
        )
        for cohort in cohorts
        for chrom in chromosomes
    )


def _predict_base_layer(
    dataset: Path,
    output_dir: Path,
    chromosomes: List[int],
    cohorts: List[str],
    cpu_count: int,
) -> None:
    """Helper function to predict base layer."""
    Parallel(n_jobs=max(1, cpu_count - 2), prefer="threads")(
        delayed(predict_window)(
            data_file=data_file,
            base_layer_model=output_dir / "base_layer_model" / "model" / f"chr{chrom}",
            output_dir=output_dir / "base_inference" / cohort / f"chr{chrom}",
            training_set=cohort == "train",
        )
        for cohort in cohorts
        for chrom in chromosomes
        for data_file in (dataset / cohort / f"chr{chrom}").rglob("data.parquet")
    )


def _save_model(
    simulated_data_path: Path,
    output_dir: Path,
    version: str,
    chromosomes: List[int],
) -> None:
    """Helper function to save the model."""
    AncestryModel(
        population_map_file=simulated_data_path / "population_map.tsv",
        model_path=output_dir.parent / version,
        sub_model_name="chr" + ".".join(map(str, chromosomes)),
        windows_info_file=output_dir / "base_layer_model" / "windows_info.tsv",
        base_model_path=output_dir / "base_layer_model" / "model",
        smooth_model_path=output_dir / "smooth_layer_model",
    ).save()


def _evaluate_results(output_dir: Path) -> None:
    """Helper function to evaluate model results."""
    # Evaluate base test
    evaluate(
        input_dir=output_dir / "base_inference" / "test",
        output_dir=output_dir / "evaluation" / "base",
        is_recomb=True,
    )
    # Evaluate final test
    evaluate(
        input_dir=output_dir / "smooth_inference_test",
        output_dir=output_dir / "evaluation" / "smooth",
    )


if __name__ == "__main__":
    training()
