"""
Misc. helper functions.
"""

from pathlib import Path
import os
import shutil
import subprocess

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


def unzip_af2_prediction(zip_file: Path, folder: Path, cleanup: bool = True) -> None:
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

    if cleanup:
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
    except FileNotFoundError:
        pass

    try:
        os.remove(unzip_dir / "cite.bibtex")
    except FileNotFoundError:
        pass

    for image_file in unzip_dir.glob("*.png"):
        try:
            os.remove(image_file)
        except FileNotFoundError:
            pass


def get_average_metric(
    results_dir: Path, metric: str = "plddt", top_only: bool = True
) -> float:
    """
    If `top_only`, get the average `metric` of the top ranked AF2 model.

    Otherwise return the average `metric` across all 5 models.
    """
    try:
        if top_only:
            result_fp = list(results_dir.glob("*scores_rank_001*.json"))[0]
            result_df = pd.read_json(result_fp)
            if metric == "pae":
                result_df["mean_residue_pae"] = result_df["pae"].apply(np.mean)
                mean_metric = np.mean(result_df["mean_residue_pae"])
            else:
                mean_metric = np.mean(result_df[metric])
        else:
            means = []
            for result_fp in results_dir.glob("*scores_rank*.json"):
                result_df = pd.read_json(result_fp)
                if metric == "pae":
                    result_df["mean_residue_pae"] = result_df["pae"].apply(np.mean)
                    means.append(np.mean(result_df["mean_residue_pae"]))
                else:
                    means.append(np.mean(result_df[metric]))
            mean_metric = np.mean(means)
    except IndexError:
        raise FileNotFoundError(
            f"This does not appear to be a directory of AlphaFold2 results. {results_dir}"
        )
    except UnicodeDecodeError as exception:
        print(f"Corrupted file for {result_fp} in {results_dir}")
        raise exception
    except ValueError as exception:
        print(f"Corrupted file for {result_fp} in {results_dir}")
        raise exception

    return mean_metric


def get_pdb_by_rank(result_dir: Path, rank: int = 1) -> Path:
    """
    Return the AF2 predicted structure by rank
    """
    pdbs = list(result_dir.glob("*_rank_001_*.pdb"))

    if len(pdbs) == 1:
        return pdbs[0]
    if len(pdbs) == 0:
        raise ValueError(
            f"Could not find any PDBs matching rank {rank} in {result_dir}"
        )

    raise ValueError(f"Found more than one PDB matching rank {rank} in {result_dir}")
