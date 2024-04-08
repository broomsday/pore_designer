"""
Define project paths.
"""

from pathlib import Path
import os


ROOT_DIR = Path(__file__).resolve().parents[1]

DATA_DIR = ROOT_DIR / "data"


def get_input_pdb_path(config: dict) -> Path:
    """
    Return the path to the input PDB used for design.
    """
    if (Path(config["directory"]) / "input" / "input.pdb").is_file():
        return Path(config["directory"]) / "input" / "input.pdb"
    elif Path(config["positive_pdb"]).is_file():
        return Path(config["positive_pdb"])

    raise FileNotFoundError("Cannot locate input PDB file")


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


def get_proteinmpnn_conda_source_path() -> Path | None:
    """
    Use an environment variable to determine where the ProteinMPNN is located.
    """
    try:
        return os.environ["PROTEINMPNN_CONDA_SOURCE"]
    except KeyError:
        print(
            "You need to set the 'PROTEINMPNN_CONDA_SOURCE' environment variable to point to the conda sourch script"
        )
        quit()


def get_proteinmpnn_env_path() -> Path | None:
    """
    Use an environment variable to determine the name of the conda environment for ProteinMPNN.
    """
    try:
        return os.environ["PROTEINMPNN_ENV"]
    except KeyError:
        print(
            "You need to set the 'PROTEINMPNN_ENV' environment variable to point to the conda environment for ProteinMPNN"
        )
        quit()


def get_alphafold_input_path(config: dict, phase: str) -> Path:
    """
    Return the path to the alphafold input.csv
    """
    return Path(config["directory"]) / "input" / f"AF2_{phase}_input.csv"


def get_alphafold_conda_source_path() -> Path | None:
    """
    Use an environment variable to determine where the alphafold conda source is located.
    """
    try:
        return os.environ["ALPHAFOLD_CONDA_SOURCE"]
    except KeyError:
        print(
            "You need to set the 'ALPHAFOLD_CONDA_SOURCE' environment variable to point to the conda sourch script for colabfold"
        )
        quit()


def get_alphafold_env_path() -> Path | None:
    """
    Use an environment variable to determine the name of the Alphafold conda environment.
    """
    try:
        return os.environ["ALPHAFOLD_ENV"]
    except KeyError:
        print(
            "You need to set the 'ALPHAFOLD_ENV' environment variable to point to the conda environment for Alphafold"
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
