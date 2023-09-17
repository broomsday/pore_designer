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
