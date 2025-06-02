"""
Utility script to generate the inputs for doing hydrophilic design with ProtienMPNN.
"""

from pathlib import Path
from shutil import copy as file_copy

import typer
import biotite.structure as bts

from designer.proteinmpnn import make_symmetry_dict, make_aa_bias_dict, make_symmetric_fixed_positions_dict, save_proteinmpnn_jsonl_dict
from designer.pdb import load_pdb


def generate_proteinmpnn_manual_run_command(input_dir: Path, pdb_name: str, num_designs: int, fixed_positions: bool, biased_aas: bool, output_dir: Path) -> str:
    """
    Generate the string that would be used to run ProteinMPNN for hydrophilic redesign.
    """
    command = f"protein_mpnn_run.py --use_soluble_model --num_seq_per_target {num_designs} --pdb_path {(input_dir / pdb_name).absolute()}"
    command += f" --tied-positions_jsonl {(input_dir / 'tied_positions.jsonl').absolute()} --out_folder {output_dir.absolute()}"

    if fixed_positions:
        command += f" --fixed-positions_jsonl {(input_dir / 'fixed_positions.jsonl').absolute()}"

    if biased_aas:
        command += f" --bias_AA_jsonl {(input_dir / 'biased_aas.jsonl').absolute()}"

    return command


def main(
        input_pdb: Path = typer.Argument(..., help="PDB that will be designed."),
        setup_dir: Path = typer.Argument(..., help="Where to copy the PDB, and save all the jsonl files."),
        fixed_positions_string: str = typer.Option(default="3,7,9,14,18", help="List of positions to fix in each chain"),
        biased_aas_string: str = typer.Option(default="V,I,L,M", help="List of amino acids to bias against"),
        bias_strength: float = typer.Option(default=-1.0, help="Strength of the bias, -ve reduces odds of design"),
        num_designs: int = typer.Option(default=10000),
    ):
    """
    Generate inputs needed to run ProteinMPNN for hydrophilic redesign.
    """
    setup_dir.mkdir(exist_ok=True, parents=True)
    input_dir = setup_dir / "inputs"
    input_dir.mkdir(exist_ok=True)
    output_dir = setup_dir / "outputs"
    output_dir.mkdir(exist_ok=True)

    structure = load_pdb(input_pdb)
    file_copy(input_pdb, input_dir / input_pdb.name)

    tied_positions_dict = make_symmetry_dict(input_pdb)
    save_proteinmpnn_jsonl_dict(tied_positions_dict, input_dir / "tied_positions.jsonl")

    if fixed_positions_string:
        chains = [chain.chain_id[0] for chain in bts.chain_iter(structure)]
        fixed_positions = [position.strip() for position in fixed_positions_string.split(",")]
        fixed_positions_dict = make_symmetric_fixed_positions_dict(input_pdb.stem, chains, fixed_positions)
        save_proteinmpnn_jsonl_dict(fixed_positions_dict, input_dir / "fixed_positions.jsonl")

    if biased_aas_string:
        biased_aas = [amino_acid.strip() for amino_acid in biased_aas_string.split(",")]
        aa_bias_dict = make_aa_bias_dict(biased_aas, bias_strength)
        save_proteinmpnn_jsonl_dict(aa_bias_dict, input_dir / "biased_aas.jsonl")
    
    print(generate_proteinmpnn_manual_run_command(input_dir, input_pdb.name, num_designs, bool(fixed_positions_string), bool(biased_aas_string), output_dir))


if __name__ == "__main__":
    typer.run(main)