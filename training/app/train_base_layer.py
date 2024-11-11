import logging
import shutil
from pathlib import Path

import pandas as pd


def get_data_from_path(path: Path, data_name: str):
    return path.parts[-1][len(data_name) :]


def train_base_layer(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate window info file
    df = pd.DataFrame(columns=["chrom", "window", "first_pos", "last_pos"])

    # Handle each chromosome directory
    for chrom_path in input_dir.glob("chr*"):
        try:
            chrom = int(get_data_from_path(chrom_path, "chr"))
        except Exception:
            logging.warning(f"{chrom_path} skipped")
            continue

        to_dir = output_dir / "model" / f"chr{chrom}"
        to_dir.mkdir(parents=True, exist_ok=True)

        # Copy generation 0 using shutil for better cross-platform support
        src_gen0 = chrom_path / "gen0"
        if src_gen0.exists():
            shutil.copytree(src_gen0, to_dir, dirs_exist_ok=True)

        # Process each window directory
        for window_path in to_dir.glob("window-*"):
            window = int(get_data_from_path(window_path, "window-"))
            df_window = pd.read_csv(window_path / "snps.tsv", sep="\t")
            new_row = {
                "chrom": chrom,
                "window": window,
                "first_pos": df_window["pos"].min(),
                "last_pos": df_window["pos"].max(),
            }
            df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

    # Final processing of the DataFrame
    df = df.sort_values(["chrom", "window"])
    df["size"] = df["last_pos"] - df["first_pos"]

    df.to_csv(output_dir / "windows_info.tsv", sep="\t", index=False)

    logging.info(f"Output dir: {output_dir}")
