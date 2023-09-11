"""
Read a 'pore designer' config file to launch of continue a pore design job.
"""


from pathlib import Path
import json

import typer
from biotite.structure.io import save_structure

from designer import pdb


app = typer.Typer()


@app.command()
def make_config(
    output_dir: Path = typer.Argument(
        ..., help="Directory to save config and job output"
    ),
    input_pdb: Path = typer.Argument(..., help="PDB to design into a pore."),
    num_mpnn: int = typer.Option(
        default=10000, help="Number of ProteinMPNN designs to make"
    ),
    num_af2: int = typer.Option(
        default=100, help="Number of design sequences to test by AlphaFold2"
    ),
    select_plddt: float = typer.Option(
        default=0.9, help="AF2 pLDDT cutoff to select final sequences"
    ),
    select_identity: float = typer.Option(
        default=0.9, help="Maximum identity between any two selected sequences"
    ),
    oligomer_rank: int = typer.Option(
        default=1,
        help="Minimum oligomer check rank to be selected for oligomer designs",
    ),
    overwrite: bool = typer.Option(default=False, help="Overwrite existing config"),
):
    """
    Based on CLI input, generate a config file for a given project.
    """
    if output_dir.is_dir():
        if not overwrite:
            raise FileExistsError("Config already  exists, use --overwrite")
        else:
            print("Directory already exists, clean manually if you want a fresh run")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "input").mkdir(exist_ok=True)
    (output_dir / "output").mkdir(exist_ok=True)
    (output_dir / "output" / "mpnn").mkdir(exist_ok=True)
    (output_dir / "output" / "af2").mkdir(exist_ok=True)
    (output_dir / "output" / "oligomer").mkdir(exist_ok=True)

    clean_input_structure = pdb.clean_pdb(input_pdb)
    clean_input_pdb = output_dir / "input" / "input.pdb"
    save_structure(clean_input_pdb, clean_input_structure)

    values = {
        "directory": str(output_dir.absolute()),
        "input_pdb": str(clean_input_pdb.absolute()),
        "num_mpnn": num_mpnn,
        "num_af2": num_af2,
        "select_plddt": select_plddt,
        "select_identity": select_identity,
        "oligomer_rank": oligomer_rank,
        "multimer": pdb.get_multimer_state(clean_input_structure),
    }

    with open(output_dir / "config.json", mode="w", encoding="utf-8") as config_file:
        json.dump(values, config_file)


@app.command()
def design_pore(
    config_file: Path = typer.Argument(..., help="Pore Design Config File")
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job.
    """


if __name__ == "__main__":
    app()
