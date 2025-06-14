"""
Unzip all AF2 prediction directories and clean them up.
"""

from pathlib import Path

import typer
from tqdm import tqdm

from designer.file_utils import unzip_af2_prediction


def main(
    in_dir: Path = typer.Argument(..., help="Directory full of AF2 .zip files"),
):
    """
    Unzip and cleanup all the .zip files.
    """
    for zip_file in tqdm(list(in_dir.glob("*result.zip")), desc="Unzipping"):
        result_dir = in_dir / zip_file.stem.replace(".result", "")
        result_dir.mkdir(exist_ok=True)
        unzip_af2_prediction(zip_file, result_dir, cleanup=True)


if __name__ == "__main__":
    typer.run(main)
