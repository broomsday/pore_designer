"""
Utility script to generate the inputs for doing hydrophilic design with ProtienMPNN.
"""

from pathlib import Path

import typer
import biotite.structure as bts
from biotite.structure.io import save_structure

from designer.proteinmpnn import make_symmetry_dict, make_aa_bias_dict, make_symmetric_fixed_positions_dict, save_proteinmpnn_jsonl_dict, rekey_proteinmpnn_dict
from designer.pdb import clean_pdb, get_sequence


def generate_proteinmpnn_manual_run_command(input_dir: Path, pdb_name: str, num_designs: int, fixed_positions: bool, biased_aas: bool, output_dir: Path) -> str:
    """
    Generate the string that would be used to run ProteinMPNN for hydrophilic redesign.
    """
    command = f"python protein_mpnn_run.py --use_soluble_model --num_seq_per_target {num_designs} --pdb_path {(input_dir / pdb_name).absolute()}"
    command += f" --tied_positions_jsonl {(input_dir / 'tied_positions.jsonl').absolute()} --out_folder {output_dir.absolute()}"

    if fixed_positions:
        command += f" --fixed_positions_jsonl {(input_dir / 'fixed_positions.jsonl').absolute()}"

    if biased_aas:
        command += f" --bias_AA_jsonl {(input_dir / 'biased_aas.jsonl').absolute()}"

    command += " --sampling_temp 0.2"

    return command


def main(
        input_pdb: Path = typer.Argument(..., help="PDB that will be designed."),
        setup_dir: Path = typer.Argument(..., help="Where to copy the PDB, and save all the jsonl files."),
        mutable_positions_string: str = typer.Option(default="3,7,10,14,18", help="List of positions to fix in each chain"),
        biased_aas_string: str = typer.Option(default="V,I,L,M,F", help="List of amino acids to bias against"),
        bias_strength: float = typer.Option(default=-1.5, help="Strength of the bias, -ve reduces odds of design"),
        num_designs: int = typer.Option(default=100),
    ):
    """
    Generate inputs needed to run ProteinMPNN for hydrophilic redesign.
    """
    setup_dir.mkdir(exist_ok=True, parents=True)
    input_dir = setup_dir / "inputs"
    input_dir.mkdir(exist_ok=True)
    output_dir = setup_dir / "outputs"
    output_dir.mkdir(exist_ok=True)

    structure = clean_pdb(input_pdb)
    save_structure(input_dir / "input.pdb", structure)

    tied_positions_dict = make_symmetry_dict(input_pdb)
    save_proteinmpnn_jsonl_dict(tied_positions_dict, input_dir / "tied_positions.jsonl")
    rekey_proteinmpnn_dict(input_dir / "tied_positions.jsonl")

    if mutable_positions_string:
        chain_lengths = [len(get_sequence(chain, mode="single")) for chain in bts.chain_iter(structure)]
        if len(set(chain_lengths)) != 1:
            raise NotImplementedError("Does not work for inputs with chains of differing lengths")

        chain_ids = [chain.chain_id[0] for chain in bts.chain_iter(structure)]
        mutable_positions = [int(position.strip()) for position in mutable_positions_string.split(",")]
        fixed_positions = [i for i in range(1, chain_lengths[0] + 1) if i not in mutable_positions]
        fixed_positions_dict = make_symmetric_fixed_positions_dict(input_pdb.stem, chain_ids, fixed_positions)
        save_proteinmpnn_jsonl_dict(fixed_positions_dict, input_dir / "fixed_positions.jsonl")
        rekey_proteinmpnn_dict(input_dir / "fixed_positions.jsonl")

    if biased_aas_string:
        biased_aas = [amino_acid.strip() for amino_acid in biased_aas_string.split(",")]
        aa_bias_dict = make_aa_bias_dict(biased_aas, bias_strength)
        save_proteinmpnn_jsonl_dict(aa_bias_dict, input_dir / "biased_aas.jsonl")


    print(generate_proteinmpnn_manual_run_command(input_dir, "input.pdb", num_designs, bool(mutable_positions_string), bool(biased_aas_string), output_dir))


if __name__ == "__main__":
    typer.run(main)