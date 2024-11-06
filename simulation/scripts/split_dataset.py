import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split


def _read_sample_map(sample_map):
    """Read and filter a sample map file.

    Args:
        sample_map (Path): Path to the sample map TSV file.

    Returns:
        pd.DataFrame: Filtered dataframe containing sample information, excluding ADMIXED
            and 'na' samples, as well as rows with null level_3_name.
    """
    df = pd.read_csv(sample_map, sep="\t")
    df["sample_id"] = df.sample_id.astype(str)
    df = df[~df.level_3_name.isin(["ADMIXED", "na"])]
    df = df[~df.level_3_name.isnull()]

    return df


def _get_stratified_sample(df, column, sample_size):
    """Get a stratified sample from the dataset based on a specified column.

    Args:
        df (pd.DataFrame): Input dataframe to sample from.
        column (str): Column name to stratify by.
        sample_size (int): Total desired sample size.

    Returns:
        pd.DataFrame: Stratified sample where each category in the specified column
            has approximately equal representation, limited by the smallest category size.
    """
    rem_size = sample_size
    dfs = []

    value_counts = df[column].value_counts(sort=True, ascending=True)
    n_categories = len(value_counts)

    for i in range(n_categories):
        count = value_counts.iloc[i]

        n_samples = int(min(count, rem_size / (n_categories - i)))
        rem_size -= n_samples

        # print(value_counts.index[i], count, n_samples)
        dfs.append(df[df[column] == value_counts.index[i]].sample(n_samples))

    return pd.concat(dfs)


def _balance_dataset(df, upper_limit=200, level_4_threshold=4):
    """Balance the dataset by removing small populations and limiting large ones.

    Args:
        df (pd.DataFrame): Input dataframe to balance.
        upper_limit (int, optional): Maximum number of samples per level 3 population.
            Defaults to 200.
        level_4_threshold (int, optional): Minimum number of samples required for a
            level 4 population. Defaults to 4.

    Returns:
        pd.DataFrame: Balanced dataframe where populations meet the size thresholds.
    """
    # trim level 4 below 4 samples
    vcs = df.meta_1.value_counts()
    sids = []

    for pop, count in vcs.items():
        if count >= level_4_threshold:
            sids += df[df.meta_1 == pop].sample_id.tolist()
        else:
            print(pop, count)
            pass  # don't include

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


def split_dataset(sample_map: str, description: str, output_dir: Path):
    """Split a sample map dataset into training and test sets.

    Args:
        sample_map (Path): Path to the input sample map TSV file.
        description (str): Description of the dataset split.
        output_dir (Path): Directory where the split datasets and metadata will be saved.

    The function:
    1. Reads and filters the sample map
    2. Balances the dataset by removing small populations and limiting large ones
    3. Performs a stratified train/test split (60/40)
    4. Saves the split datasets and metadata
    """
    # read sample map
    sample_map_df = _read_sample_map(sample_map)
    sample_map_df = sample_map_df[sample_map_df.level_3_code != "na"]

    print(f"Number of populations: {len(sample_map_df.level_3_name.drop_duplicates())}")
    print(f"Number of samples: {len(sample_map_df.sample_id)}")
    print(f"Number of samples (after level 1 filter): {len(sample_map_df.sample_id)}")

    # balance dataset by trimming and removing populations
    df = _balance_dataset(sample_map_df, upper_limit=500)
    # perform train/test split
    df_train, df_test = train_test_split(
        df,
        test_size=0.4,
        shuffle=True,
        stratify=df[["level_3_name", "meta_1"]],
    )

    print(f"Train Samples: {len(df_train)}")
    print(f"Test Samples : {len(df_test)}")

    assert len(df_train) + len(df_test) == len(sample_map_df)

    # make sure test/train split doesn't overlap
    assert (
        np.intersect1d(df_test.sample_id.to_numpy(), df_train.sample_id.to_numpy()).size
        == 0
    )

    df_train.to_csv(output_dir / "train/pheno/sample_map.tsv", sep="\t", index=False)
    df_test.to_csv(output_dir / "test/pheno/sample_map.tsv", sep="\t", index=False)

    with open(output_dir / "meta.json", "w") as f:
        json.dump(
            {
                "sample_map": sample_map,
                "description": description,
            },
            f,
        )

    df = (
        sample_map_df[["level_3_code", "level_1_name", "level_2_name", "level_3_name"]]
        .drop_duplicates()
        .sort_values(["level_1_name", "level_2_name", "level_3_name"])
        .reset_index()
        .drop("index", axis=1)
        .reset_index()
        .rename(columns={"index": "id"})
    )
    df = df[["id", "level_3_code", "level_1_name", "level_2_name", "level_3_name"]]
    df.to_csv(
        output_dir / "population_map.tsv",
        sep="\t",
        index=False,
    )

    return {
        "train": output_dir / "train/pheno/sample_map.tsv", 
        "test": output_dir / "test/pheno/sample_map.tsv",
    }
