"""
Misc. helper functions.
"""


from pathlib import Path
import os
import shutil
import subprocess
import json

import numpy as np
import pandas as pd


EXPECTED_NUM_PREDICTIONS = 5


def keep_files_by_suffix(folder: Path, suffixes: list[str]) -> None:
    """
    Delete all files in `folder` that do not end with a suffix in `suffixes`
    """
    for file_path in folder.glob("*.*"):
        if file_path.suffix not in suffixes:
            os.remove(file_path)


def unzip_af2_prediction(zip_file: Path, folder: Path) -> None:
    """
    Unzip a file into a new directory
    """
    folder.mkdir(exist_ok=True, parents=True)

    # don't unzip if we already have the individual files
    if len(list(folder.glob("*.json"))) >= EXPECTED_NUM_PREDICTIONS:
        return None

    # copy over the archive
    file_name = f"{folder.parts[-1]}.result.zip"
    copy_archive = folder / file_name
    if not copy_archive.is_file():
        original_archive = zip_file
        shutil.copy(original_archive, copy_archive)

    # unzip
    cwd = os.getcwd()
    os.chdir(folder)
    subprocess.run(["unzip", file_name], stdout=subprocess.DEVNULL, check=False)
    os.chdir(cwd)

    # cleanup
    clean_af2_prediction_files(folder)


def clean_af2_prediction_files(unzip_dir: Path) -> None:
    """
    Remove the archive, citation, and images.
    Remove all but the top-ranked PDB model.
    Reduce the .json files to only the pLDDT for space.
    """
    file_stem = f"{unzip_dir.parts[-1]}"

    # remove unnecessary files
    try:
        os.remove(unzip_dir / f"{file_stem}.result.zip")
        os.remove(unzip_dir / "cite.bibtex")
    except FileNotFoundError:
        pass  # zip only present if we copied it over, bibtex sometimes missing
    os.remove(unzip_dir / f"{file_stem}_predicted_aligned_error_v1.json")

    for image_file in unzip_dir.glob("*.png"):
        os.remove(image_file)

    # reduce .json files to just the pLDDT entries
    for scores_file in unzip_dir.glob(f"{file_stem}*.json"):
        with open(scores_file, mode="r", encoding="utf-8") as fi:
            scores = json.load(fi)
        plddts = {"plddt": scores["plddt"]}

        with open(scores_file, mode="w", encoding="utf-8") as fo:
            json.dump(plddts, fo)


def get_average_plddt(results_dir: Path, top_only: bool = True) -> float:
    """
    If `top_only`, get the average pLDDT of the top ranked AF2 model.

    Otherwise return the average pLDDT across all 5 models.
    """
    try:
        if top_only:
            result_fp = list(results_dir.glob("*scores_rank_001*.json"))[0]
            result_df = pd.read_json(result_fp)
            mean_plddt = np.mean(result_df["plddt"])
        else:
            means = []
            for result_fp in results_dir.glob("*scores_rank*.json"):
                result_df = pd.read_json(result_fp)
                means.append(np.mean(result_df["plddt"]))
            mean_plddt = np.mean(means)
    except IndexError:
        raise FileNotFoundError(
            f"This does not appear to be a directory of AlphaFold2 results. {results_dir}"
        )

    return mean_plddt
