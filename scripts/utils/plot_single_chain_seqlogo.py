"""
Plot a sequence logo of the first chain, for all sequences in an MPNN fasta, or all PDBs in a directory.
"""

from pathlib import Path

import typer

from designer.plotting import seqlogo


def main(
    in_mpnn_fasta: Path = typer.Argument(..., help="ProteinMPNN result FASTA."),
    out_logo_dir: Path = typer.Argument(..., help="Where to save the sequence logos."),
    top_num: int = typer.Option(default=100, help="How many top designs by MPNN score to use for the second sequence logo."),
):
    """
    Plot a sequence logo of all designs and the `top_num` designs
    """


if __name__ == "__main__":
    typer.run(main)