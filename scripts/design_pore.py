"""
Read a 'pore designer' config file to launch of continue a pore design job.
"""


from pathlib import Path
import shutil
import json

import typer
from biotite.structure.io import save_structure

from designer import pdb, proteinmpnn, alphafold, utils, paths


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
    temperature_mpnn: float = typer.Option(
        default=0.1, help="MPNN sampling temperature."
    ),
    num_af2: int = typer.Option(
        default=100, help="Number of design sequences to test by AlphaFold2"
    ),
    recycle_af2: int = typer.Option(default=3, help="Number of recycles to use in AF2"),
    select_plddt: float = typer.Option(
        default=0.9, help="AF2 pLDDT cutoff to select final sequences"
    ),
    select_rmsd: float = typer.Option(default=1.0, help="Max RMSD for selection"),
    select_identity: float = typer.Option(
        default=0.9, help="Maximum identity between any two selected sequences"
    ),
    select_oligomer_rank: int = typer.Option(
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

    clean_input_structure = pdb.clean_pdb(input_pdb)
    clean_input_pdb = output_dir / "input" / "input.pdb"
    save_structure(clean_input_pdb, clean_input_structure)

    values = {
        "directory": str(output_dir.absolute()),
        "input_pdb": str(clean_input_pdb.absolute()),
        "num_mpnn": num_mpnn,
        "temperature_mpnn": temperature_mpnn,
        "num_af2": num_af2,
        "recycle_af2": recycle_af2,
        "select_plddt": select_plddt,
        "select_rmsd": select_rmsd,
        "select_identity": select_identity,
        "select_oligomer_rank": select_oligomer_rank,
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
        proteinmpnn.rename_existing_results(config)
        if proteinmpnn.get_num_to_design(config) > 0:
            proteinmpnn_script = proteinmpnn.make_shell_script(config)
            proteinmpnn.run_proteinmpnn(proteinmpnn_script)

        proteinmpnn_seqs = proteinmpnn.select_top_sequences(config)
        proteinmpnn_seqs_dict = [seq._asdict() for seq in proteinmpnn_seqs]
        proteinmpnn.save_top_sequences(config, proteinmpnn_seqs_dict)
    proteinmpnn_seqs = proteinmpnn.load_top_sequences(config)

    # setup AF2 input file
    completed_af2_ids = alphafold.get_completed_ids(config, "design")
    if len(completed_af2_ids) < config["num_af2"]:
        af2_input_df = alphafold.make_af2_design_input(
            proteinmpnn_seqs, completed_af2_ids
        )
        af2_input_df.to_csv(paths.get_alphafold_input_path(config, "design"))

        # run AF2 batch and process the output
        alphafold_script = alphafold.make_shell_script(config, "design")
        alphafold.run_af2_batch(alphafold_script)

    if not alphafold.have_selected(config, "design"):
        # unzip everything, parse pLDDT, compute RMSD to the templates and filter
        selected = alphafold.select_top(
            config,
            plddt=config["select_plddt"],
            rmsd=config["select_rmsd"],
            oligomer=None,
            identity=config["select_identity"],
        )
        utils.save_selected(config, "design", selected)
    selected = utils.load_selected(config, "design")

    if config["multimer"] > 1:
        if not alphafold.have_selected(config, "oligomer"):
            # generate the oligomer check inputs
            af2_input_df = alphafold.make_af2_oligomer_input(selected)
            af2_input_df.to_csv(paths.get_alphafold_input_path(config, "oligomer"))

            # run the AF2 oligomer check
            alphafold_script = alphafold.make_shell_script(config, "oligomer")
            alphafold.run_af2_batch(alphafold_script)

            # pick the winners of the oligomer check
            selected = alphafold.select_top(
                config,
                plddt=config["select_plddt"],
                rmsd=config["select_rmsd"],
                oligomer=config["select_oligomer_rank"],
                identity=config["select_identity"],
            )
            utils.save_selected(config, "oligomer", selected)
        selected = utils.load_selected(config, "oligomer")

    # report the winners and copy files
    utils.report_selected(config, selected)

    print(config)


if __name__ == "__main__":
    app()
