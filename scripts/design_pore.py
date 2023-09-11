"""
Read a 'pore designer' config file to launch of continue a pore design job.
"""


from pathlib import Path

import typer


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
        defalt=100, help="Number of design sequences to test by AlphaFold2"
    ),
    select_plddt: float = typer.Option(
        default=0.9, help="AF2 pLDDT cutoff to select final sequences"
    ),
    select_identity: float = typer.Option(
        default=0.9, help="Maximum identity between any two selected sequences"
    ),
):
    """
    Based on CLI input, generate a config file for a given project.
    """


@app.command()
def design_pore(
    config_file: Path = typer.Argument(..., help="Pore Design Config File")
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job.
    """
