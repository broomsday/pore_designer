"""
Define project paths.
"""

from pathlib import Path
import os


ROOT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT_DIR / "data"


def get_proteinmpnn_path() -> Path | None:
    """
    Use an environment variable to determine where ProteinMPNN is installed.
    """
    try:
        return os.environ["PROTEINMPNN"]
    except KeyError:
        print("You need to install 'proteinmpnn'")
        print(
            "Then set the 'PROTEINMPNN' environment variable to point to it's location"
        )
        quit()


def get_conda_source_path() -> Path | None:
    """
    Use an environment variable to determine where the conda source shell script is located.
    """
    try:
        return os.environ["CONDA_SOURCE"]
    except KeyError:
        print(
            "You need to set the 'CONDA_SOURCE' environment variable to point to conda's sourch.sh script"
        )
        quit()


def get_input_pdb_path(config: dict) -> Path:
    """
    Return the path to the input PDB used for design.
    """
    return Path(config["directory"]) / "input" / "input.pdb"


def get_alphafold_input_path(config: dict, phase: str) -> Path:
    """
    Return the path to the alphafold input.csv
    """
    return Path(config["directory"]) / "input" / f"AF2_{phase}_input.csv"


def get_alphafold_source_path() -> Path | None:
    """
    Use an environment variable to determine where the alphafold conda source is located.
    """
    try:
        return os.environ["ALPHAFOLD_ENV"]
    except KeyError:
        print(
            "You need to set the 'ALPHAFOLD_ENV' environment variable to point to the conda sourch script for colabfold"
        )
        quit()


def get_alphafold_results_path(config: dict, phase: str) -> Path:
    """
    The path used for saving the AlphaFold results.
    """
    return Path(config["directory"]) / f"top_alphafold_{phase}_results.json"


def get_alphafold_selected_path(config: dict, phase: str) -> Path:
    """
    The path used for saving the AlphaFold selected sequences.
    """
    return Path(config["directory"]) / f"top_alphafold_{phase}_selected.json"
