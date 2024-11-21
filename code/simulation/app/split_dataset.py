import json
import logging
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

logging.basicConfig(level=logging.INFO)


def _read_sample_map(sample_map: str) -> pd.DataFrame:
    """Read and filter a sample map file.

    Args:
        sample_map (str): Path to the sample map TSV file.

    Returns:
        pd.DataFrame: Filtered dataframe containing:
            - Excludes samples with level_3_name of "ADMIXED" or "na"
            - Excludes rows with null level_3_name
            - Contains columns: sample_id (str), level_3_name, meta_1, etc.

    Raises:
        Exception: If there's an error reading or processing the sample map file.
    """
    try:
        df = pd.read_csv(sample_map, sep="\t", dtype={"sample_id": str})
        return df[
            (~df.level_3_name.isin(["ADMIXED", "na"])) & (df.level_3_name.notna())
        ]
    except Exception as e:
        logging.error(f"Error reading sample map: {e}")
        raise


def _get_stratified_sample(
    df: pd.DataFrame, column: str, sample_size: int
) -> pd.DataFrame:
    """Get a stratified sample from the dataset based on a specified column.

    Args:
        df (pd.DataFrame): Input dataframe containing samples to select from.
        column (str): Column name to use for stratification (e.g., "meta_1").
        sample_size (int): Target total number of samples to return.

    Returns:
        pd.DataFrame: Stratified sample with approximately equal representation per category,
            limited by the smallest category size. Contains same columns as input df.
    """
    value_counts = df[column].value_counts(sort=True, ascending=True)
    samples = []
    remaining_size = sample_size
    n_categories = len(value_counts)

    for i, (category, count) in enumerate(value_counts.items()):
        n_samples = min(count, remaining_size // (n_categories - i))
        remaining_size -= n_samples
        samples.append(df[df[column] == category].sample(n_samples))

    return pd.concat(samples, ignore_index=True)


def _balance_dataset(
    df: pd.DataFrame, upper_limit: int = 200, level_4_threshold: int = 4
) -> pd.DataFrame:
    """Balance the dataset by removing small populations and limiting large ones.

    Args:
        df (pd.DataFrame): Input dataframe containing sample information.
        upper_limit (int, optional): Maximum samples allowed per level 3 population.
            Populations exceeding this will be downsampled. Defaults to 200.
        level_4_threshold (int, optional): Minimum samples required for a level 4
            population to be included. Defaults to 4.

    Returns:
        pd.DataFrame: Balanced dataframe where:
            - Level 4 populations have at least level_4_threshold samples
            - Level 3 populations have at most upper_limit samples
            - Small populations (<4 samples) are removed
    """
    # trim level 4 below 4 samples
    vcs = df.meta_1.value_counts()
    sids = []

    for pop, count in vcs.items():
        if count >= level_4_threshold:
            sids += df[df.meta_1 == pop].sample_id.tolist()
        else:
            logging.error(f"{pop} has {count} samples (<{level_4_threshold}), skipping")

    df = df[df.sample_id.isin(sids)]

    # trim populations to size of 10-100 samples
    vcs = df.level_3_name.value_counts()
    sids = []

    for pop, count in vcs.items():
        if count > upper_limit:
            sids += _get_stratified_sample(
                df[df.level_3_name == pop],
                "meta_1",
                upper_limit,
            ).sample_id.tolist()
        elif count >= 4:
            sids += df[df.level_3_name == pop].sample_id.tolist()
        else:
            pass  # don't include

    df = df[df.sample_id.isin(sids)]

    return df


def split_dataset(
    sample_map: str, description: str, temp_path: Path, output_path: Path
) -> Dict[str, Path]:
    """Split a sample map dataset into training and test sets.

    This function processes a sample map file, balances the dataset, and creates
    train/test splits while maintaining population stratification.

    Args:
        sample_map (str): Path to the input sample map TSV file.
        description (str): Version/Description of the dataset split.
        output_path (Path): Directory where the split datasets and metadata will be saved.

    Returns:
        Dict[str, Path]: Dictionary containing paths to the output files for each split
            ('train' and 'test').

    Raises:
        ValueError: If the split validation fails (size mismatch or overlapping samples).
        Exception: If there's an error during processing.

    The function performs the following steps:
        1. Reads and filters the sample map
        2. Balances the dataset by removing small populations and limiting large ones
        3. Performs a stratified train/test split (60/40)
        4. Saves the split datasets and metadata
    """
    try:
        # Read and filter data
        sample_map_df = _read_sample_map(sample_map)
        sample_map_df = sample_map_df.loc[sample_map_df.level_3_code != "na"]

        # Log initial statistics
        logging.info("Dataset statistics:")
        logging.info(f"- Populations: {sample_map_df.level_3_name.nunique()}")
        logging.info(f"- Total samples: {len(sample_map_df)}")

        # Balance dataset
        df = _balance_dataset(sample_map_df, upper_limit=500)

        # Perform split
        train_df, test_df = train_test_split(
            df,
            test_size=0.4,
            shuffle=True,
            stratify=df[["level_3_name", "meta_1"]],
        )

        # Validate split
        if not (len(train_df) + len(test_df) == len(df)):
            raise ValueError("Split validation failed: size mismatch")
        if np.intersect1d(test_df.sample_id, train_df.sample_id).size != 0:
            raise ValueError("Split validation failed: overlapping samples")

        # Save splits
        result = {}
        for cohort, cohort_df in [("train", train_df), ("test", test_df)]:
            logging.info(f"{cohort.title()} samples: {len(cohort_df)}")
            path = temp_path / cohort / "pheno"
            path.mkdir(parents=True, exist_ok=True)
            sample_map_path = path / "sample_map.tsv"
            cohort_df.to_csv(sample_map_path, sep="\t", index=False)
            result[cohort] = sample_map_path

        # Save metadata
        _save_metadata(output_path, sample_map, description, sample_map_df)

        return result

    except Exception as e:
        logging.error(f"Error in split_dataset: {e}")
        raise


def _save_metadata(
    output_path: Path, sample_map: Path, description: str, sample_map_df: pd.DataFrame
) -> None:
    """Save metadata files for the dataset split.

    Args:
        output_path (Path): Directory where metadata files will be saved.
        sample_map (Path): Original sample map file path.
        description (str): Description of the dataset split.
        sample_map_df (pd.DataFrame): The processed sample map dataframe.

    The function saves two files:
        1. meta.json: Contains basic metadata about the split
        2. population_map.tsv: Contains detailed population information
    """

    logging.info(f"Saving meta info: {output_path / 'meta.json'}")
    with open(output_path / "meta.json", "w") as f:
        data = {
            "sample_map": str(sample_map),
            "description": description,
        }
        json.dump(data, f)

    logging.info(f"Saving population map: {output_path / 'population_map.tsv'}")
    (
        sample_map_df[["level_3_code", "level_1_name", "level_2_name", "level_3_name"]]
        .drop_duplicates()
        .sort_values(["level_1_name", "level_2_name", "level_3_name"])
        .reset_index(drop=True)
        .reset_index()
        .rename(columns={"index": "id"})[
            ["id", "level_3_code", "level_1_name", "level_2_name", "level_3_name"]
        ]
        .to_csv(output_path / "population_map.tsv", sep="\t", index=False)
    )
