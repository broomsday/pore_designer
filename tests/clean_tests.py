"""
Clean the test directories so they can be reused.
"""

import os
import shutil

from designer.paths import ROOT_DIR


TESTS = [
    "metric",
    "positive",
    "negative",
    "mutation",
]

DIRECTORIES = ["output"]
PROTECTED_PATHS = ["config.json", "input"]


def clean_test_dirs():
    """
    Clean the test directories of output so they are ready for a fresh run.
    """
    for test in TESTS:
        top_dir = ROOT_DIR / "tests" / f"{test}_test"

        # perform general cleaning
        for directory in DIRECTORIES:
            clean_dir = top_dir / directory
            for file_path in clean_dir.glob("*"):
                if file_path.name not in PROTECTED_PATHS:
                    if file_path.is_file():
                        os.remove(file_path)
                    elif file_path.is_dir():
                        shutil.rmtree(file_path)

        # specifically clean up the input subdirectory of output
        out_input_dir = top_dir / "output" / "input"
        for file_path in out_input_dir.glob("**/*.*"):
            if (file_path.suffix != ".pdb") and (file_path.stem != "input"):
                os.remove(file_path)


if __name__ == "__main__":
    clean_test_dirs()
