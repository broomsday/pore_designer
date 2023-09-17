"""
Read a 'pore designer' config file to launch of continue a pore design job.
"""


from pathlib import Path
import shutil
import json

import typer
from biotite.structure.io import save_structure

from designer import pdb, proteinmpnn


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
    symmetry_dict: Path = typer.Option(
        default=None,
        help="If using a monomer from 'monomerizer' provide the symmetry dict output",
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

    if symmetry_dict is not None:
        values["symmetry_dict"] = str(symmetry_dict.absolute())
    else:
        values["symmetry_dict"] = None

    with open(output_dir / "config.json", mode="w", encoding="utf-8") as config_file:
        json.dump(values, config_file)


@app.command()
def design_pore(
    config_path: Path = typer.Argument(..., help="Pore Design config file")
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job.
    """
    with open(config_path, mode="r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    # setup the symmetry dict for proteinmpnn unless we already have it
    symmetry_path = Path(config["input_pdb"]).with_name("input.symm.jsonl")
    if config["symmetry_dict"] is None:
        symmetry_dict = proteinmpnn.make_symmetry_dict(Path(config["input_pdb"]))
        proteinmpnn.save_symmetry_dict(symmetry_dict, symmetry_path)
        config["symmetry_dict"] = symmetry_path
    else:
        shutil.copy(config["symmetry_dict"], symmetry_path)
        proteinmpnn.rekey_symmetry_dict(symmetry_path)
    config["symmetry_dict"] = str(symmetry_path)

    # run proteinmpnn and get the best sequences
    if not proteinmpnn.have_top_results(config):
        if proteinmpnn.get_num_to_design(config) > 0:
            proteinmpnn_script = proteinmpnn.make_shell_script(config)
            proteinmpnn.run_proteinmpnn(proteinmpnn_script)

        proteinmpnn_seqs = proteinmpnn.select_top_sequences(config)
        proteinmpnn_seqs_dict = [seq._asdict() for seq in proteinmpnn_seqs]
        proteinmpnn.save_top_sequences(config, proteinmpnn_seqs_dict)

    proteinmpnn_seqs = proteinmpnn.load_top_sequences(config)

    print(proteinmpnn_seqs)

    # TODO: run AF2 on the top ProteinMPNN designs
    # TODO: pick the winners

    # TODO: if a multimer, run the AF2 multimer check
    # TODO: pick the winners

    # TODO: report the winners

    print(config)


if __name__ == "__main__":
    app()
