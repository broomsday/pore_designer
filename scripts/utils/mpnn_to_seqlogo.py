"""
Plot a sequence logo of the first chain, for all sequences in an MPNN fasta.
"""

from pathlib import Path

import typer

from designer.proteinmpnn import get_sequences, get_top_mpnn_sequences_by_score
from designer.plotting import seqlogo


def main(
    in_mpnn_fasta: Path = typer.Argument(..., help="ProteinMPNN result FASTA."),
    out_logo_dir: Path = typer.Argument(..., help="Where to save the sequence logos."),
    top_num: int = typer.Option(default=100, help="How many top designs by MPNN score to use for the second sequence logo."),
):
    """
    Plot a sequence logo of all designs and the `top_num` designs
    """
    out_logo_dir.mkdir(exist_ok=True, parents=True)

    mpnn_sequences = get_sequences(in_mpnn_fasta)
    seqlogo([mpnn_sequence.sequence[0] for mpnn_sequence in mpnn_sequences], out_logo_dir / "seq_logo_all.png")

    top_mpnn_sequences = get_top_mpnn_sequences_by_score(mpnn_sequences, top_num)
    seqlogo([mpnn_sequence.sequence[0] for mpnn_sequence in top_mpnn_sequences], out_logo_dir / f"seq_logo_top_{top_num}.png")
    

if __name__ == "__main__":
    typer.run(main)