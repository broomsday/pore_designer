"""
Utility script to generate the inputs for doing hydrophilic design with ProtienMPNN.
"""

from pathlib import Path
from shutil import copy as file_copy

import typer
import biotite.structure as bts

from designer.proteinmpnn import make_symmetry_dict, make_aa_bias_dict, make_symmetric_fixed_positions_dict, save_proteinmpnn_jsonl_dict
from designer.pdb import load_pdb


def main(
        input_pdb: Path = typer.Argument(..., help="PDB that will be designed."),
        output_dir: Path = typer.Argument(..., help="Where to copy the PDB, and save all the jsonl files."),
        fixed_positions_string: str = typer.Option(default="3,7,9,14,18", help="List of positions to fix in each chain"),
        biased_aas_string: str = typer.Option(default="V,I,L,M", help="List of amino acids to bias against"),
        bias_strength: float = typer.Option(default=-1.0, help="Strength of the bias, -ve reduces odds of design"),
    ):
    """
    Generate inputs needed to run ProteinMPNN for hydrophilic redesign.
    """
    output_dir.mkdir(exist_ok=True, parents=True)

    structure = load_pdb(input_pdb)
    file_copy(input_pdb, output_dir / input_pdb.name)

    tied_positions_dict = make_symmetry_dict(input_pdb)
    save_proteinmpnn_jsonl_dict(tied_positions_dict, output_dir / "tied_positions.jsonl")

    if fixed_positions_string:
        chains = [chain.chain_id[0] for chain in bts.chain_iter(structure)]
        fixed_positions = [position.strip() for position in fixed_positions_string.split(",")]
        fixed_positions_dict = make_symmetric_fixed_positions_dict(input_pdb.stem, chains, fixed_positions)
        save_proteinmpnn_jsonl_dict(fixed_positions_dict, output_dir / "fixed_positions.jsonl")

    if biased_aas_string:
        biased_aas = [amino_acid.strip() for amino_acid in biased_aas_string.split(",")]
        aa_bias_dict = make_aa_bias_dict(biased_aas, bias_strength)
        save_proteinmpnn_jsonl_dict(aa_bias_dict, output_dir / "biased_aas.jsonl")
    

if __name__ == "__main__":
    typer.run(main)