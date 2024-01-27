"""
Module to hold top-level design functions
"""


from pathlib import Path
import shutil
import json

import typer
import pandas as pd

from designer import proteinmpnn, alphafold, paths


def design_pore_positive(
    config_path: Path = typer.Argument(
        ..., help="Pore Design config file for positive design"
    )
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job for positive design.
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
        proteinmpnn.rekey_proteinmpnn_dict(symmetry_path)
    config["symmetry_dict"] = str(symmetry_path)

    # setup the fixed dict for proteinmpnn if it was provided
    if config["fixed_dict"] is not None:
        fixed_path = Path(config["input_pdb"]).with_name("input.fixed.jsonl")
        shutil.copy(config["fixed_dict"], fixed_path)
        proteinmpnn.rekey_proteinmpnn_dict(fixed_path)
        config["fixed_dict"] = str(fixed_path)

    # run proteinmpnn and get the best sequences
    if not proteinmpnn.have_top_results(config):
        proteinmpnn.rename_existing_results(config)
        if proteinmpnn.get_num_to_design(config) > 0:
            proteinmpnn_script = proteinmpnn.make_shell_script(config)
            proteinmpnn.run_proteinmpnn(proteinmpnn_script)

        proteinmpnn_seqs = proteinmpnn.select_top_sequences(config)
        proteinmpnn.save_top_sequences(config, proteinmpnn_seqs)
    proteinmpnn_seqs = proteinmpnn.load_top_sequences(config)

    # setup AF2 input file
    completed_af2_ids = alphafold.get_completed_ids(config, "design")
    if len(completed_af2_ids) < config["num_af2"]:
        af2_input_df = alphafold.make_af2_design_input(
            proteinmpnn_seqs, completed_af2_ids
        )
        af2_input_df.to_csv(paths.get_alphafold_input_path(config, "design"))

        # run AF2 batch
        alphafold_script = alphafold.make_shell_script(config, "design")
        alphafold.run_af2_batch(alphafold_script)

    # analyze alphafold results
    if not alphafold.have_alphafold(config, "design", "results"):
        design_results = alphafold.compile_alphafold_design_results(config)
        alphafold.save_alphafold(config, "design", "results", design_results)
    design_results = alphafold.load_alphafold(config, "design", "results")

    # select top sequences before any oligomer checking
    if not alphafold.have_alphafold(config, "design", "selected"):
        design_selected = alphafold.select_top(
            design_results,
            top_plddt=config["top_plddt"],
            mean_plddt=config["mean_plddt"],
            top_rmsd=config["top_rmsd"],
            mean_rmsd=config["mean_rmsd"],
            oligomer=None,
            max_identity=config["select_identity"],
        )
        alphafold.save_alphafold(config, "design", "selected", design_selected)
    design_selected = alphafold.load_alphafold(config, "design", "selected")
    print(f"Selected: {len(design_selected)} sequences from initial design")

    # if we're just doing a monomer, report now
    if config["multimer"] == 1:
        alphafold.report_selected(config, design_selected)
        quit()

    # perform the alphafold oligomer check for multimers
    completed_af2_ids = alphafold.get_completed_ids(config, "oligomer")
    if len(completed_af2_ids) < (
        len(design_selected) * (alphafold.LOWER_OLIGOMERS + alphafold.HIGHER_OLIGOMERS)
    ):
        af2_input_df = alphafold.make_af2_oligomer_input(
            design_selected, completed_af2_ids
        )
        af2_input_df.to_csv(paths.get_alphafold_input_path(config, "oligomer"))

        # run AF2 batch
        alphafold_script = alphafold.make_shell_script(config, "oligomer")
        alphafold.run_af2_batch(alphafold_script)

    # analyze alphafold results
    if not alphafold.have_alphafold(config, "oligomer", "results"):
        oligomer_results = alphafold.compile_alphafold_oligomer_results(
            config, design_selected
        )
        alphafold.save_alphafold(config, "oligomer", "results", oligomer_results)
    oligomer_results = alphafold.load_alphafold(config, "oligomer", "results")

    # select top sequences
    if not alphafold.have_alphafold(config, "oligomer", "selected"):
        oligomer_selected = alphafold.select_top(
            oligomer_results,
            top_plddt=config["top_plddt"],
            mean_plddt=config["mean_plddt"],
            top_rmsd=config["top_rmsd"],
            mean_rmsd=config["mean_rmsd"],
            oligomer=config["select_oligomer_rank"],
            max_identity=config["select_identity"],
        )
        alphafold.save_alphafold(config, "oligomer", "selected", oligomer_selected)
    oligomer_selected = alphafold.load_alphafold(config, "oligomer", "selected")
    print(f"Selected: {len(oligomer_selected)} sequences from oligomer check")

    alphafold.report_selected(config, oligomer_selected)


def design_pore_negative(
    config_path: Path = typer.Argument(
        ..., help="Pore Design config file for negative design"
    )
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job for negative design.
    """
    with open(config_path, mode="r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    # TODO: run ProteinMPNN on positive and negative

    # TODO: produce positive->negative distance metrics

    # TODO: choose best designs based on metrics

    # TODO: run AF2

    # TODO: report best designs

    print(config)


def produce_multimer_metrics(
    config_path: Path = typer.Argument(
        ..., help="Pore Design config file for metric production"
    )
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job for metric production.
    """
    with open(config_path, mode="r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    input_data = pd.read_csv(config["input_csv"], index_col="pdb")
    print(input_data)

    # TODO: evaluate metrics for all sequences at different multimer values

    # TODO: add mean PAE metric or RMS-PAE

    print(config)
