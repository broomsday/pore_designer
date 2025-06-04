"""
Generate an AF2 input.csv from the top N scoring sequences from an MPNN run.
"""

from pathlib import Path

import typer

from designer.proteinmpnn import get_sequences, get_top_mpnn_sequences_by_score
from designer.alphafold import make_af2_design_input

def main(
    in_mpnn: Path = typer.Argument(..., help="MPNN FASTA output file."),
    out_csv: Path = typer.Argument(..., help="Save location for the AF2 input.csv"),
    name: str = typer.Argument(..., help="AF2 input name for id column."),
    num_designs: int = typer.Option(default=100, help="How many top sequences to use for input"),
):
    """
    Read in the MPNN sequences and scores and pick the top `num_designs` by their MPNN score and build an AF2 input.csv of those.
    """
    mpnn_seqs = get_sequences(in_mpnn)
    top_mpnn_seqs = get_top_mpnn_sequences_by_score(mpnn_seqs, num_designs)

    af2_input = make_af2_design_input(top_mpnn_seqs, [])
    af2_input["id"] = af2_input.index.to_list()
    af2_input["id"] = af2_input["id"].apply(lambda id: id.replace("None", name))
    af2_input = af2_input.set_index("id")

    af2_input.to_csv(out_csv)
    print(af2_input)


if __name__ == "__main__":
    typer.run(main)