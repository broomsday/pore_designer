"""
Module to hold top-level design functions
"""

from pathlib import Path
import shutil
import json

import typer
import pandas as pd

from designer import proteinmpnn, alphafold, paths, plotting


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

    # run ProteinMPNN on positive
    positive_pdb = Path(config["positive_pdb"])
    symmetry_dict = proteinmpnn.make_symmetry_dict(positive_pdb)
    symmetry_path = Path(positive_pdb).with_name(f"{positive_pdb.stem}.symm.jsonl")
    proteinmpnn.save_symmetry_dict(symmetry_dict, symmetry_path)
    config["symmetry_dict"] = symmetry_path

    if proteinmpnn.get_num_to_design(config, positive_pdb) > 0:
        proteinmpnn.rename_existing_results(config, positive_pdb)
        proteinmpnn_script = proteinmpnn.make_shell_script(config, positive_pdb)
        proteinmpnn.run_proteinmpnn(proteinmpnn_script)

    # run ProteinMPNN on negatives
    for negative_pdb in Path(config["negative_pdbs"]).glob("*.pdb"):
        symmetry_dict = proteinmpnn.make_symmetry_dict(negative_pdb)
        symmetry_path = Path(negative_pdb).with_name(f"{negative_pdb.stem}.symm.jsonl")
        proteinmpnn.save_symmetry_dict(symmetry_dict, symmetry_path)
        config["symmetry_dict"] = symmetry_path

        if proteinmpnn.get_num_to_design(config, negative_pdb) > 0:
            proteinmpnn.rename_existing_results(config, negative_pdb)
            proteinmpnn_script = proteinmpnn.make_shell_script(config, negative_pdb)
            proteinmpnn.run_proteinmpnn(proteinmpnn_script)

    print(config)

    # TODO: summarize positive and negative into distributions

    # TODO: compute difference distributions for positive to each negative

    # TODO: score each positive sequence by similarity to each difference distribution

    # TODO: compute the overall similarity and min-similarity metric

    # TODO: choose best designs based on metrics (half by highest overall similarity, half by higher min similarity)

    # TODO: run AF2 as usual on WT and Oligomers for the sequences selected above

    # TODO: report best designs


def metric_test_to_select_seq(
    pdb: str, sequence: str, oligomer: int
) -> alphafold.SelectSeq:
    """
    Generate a SelectSeq object from a pdb_id, sequence, and oligomer count.
    """
    return alphafold.SelectSeq(
        id=pdb,
        sequence=[sequence] * oligomer,
        merged_sequence=":".join([sequence] * oligomer),
        unique_chains=len(sequence.split(":")),
        score=None,
        recovery=None,
        source="wt",
        mutation=None,
        frequency=None,
        selection=None,
        top_plddt=None,
        mean_plddt=None,
        top_rmsd=None,
        mean_rmsd=None,
        top_oligomer=None,
        designed_oligomer_rank=None,
    )


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

    # build the input data into SelectSeq objects
    input_data = pd.read_csv(config["input_csv"], index_col="pdb")
    metric_seqs = [
        metric_test_to_select_seq(
            input_data.index[idx],
            input_data.iloc[idx]["sequence"],
            input_data.iloc[idx]["oligomer"],
        )
        for idx in range(len(input_data.index))
    ]

    # run AF2 on the WT oligomer
    completed_af2_ids = alphafold.get_completed_ids(config, "design")
    if len(completed_af2_ids) < len(metric_seqs):
        af2_input_df = alphafold.make_af2_wt_metric_input(
            metric_seqs, completed_af2_ids
        )
        af2_input_df.to_csv(paths.get_alphafold_input_path(config, "design"))

        # run AF2 batch
        alphafold_script = alphafold.make_shell_script(config, "design")
        alphafold.run_af2_batch(alphafold_script)

    # run AF2 on the expanded oligomer cases
    completed_af2_ids = alphafold.get_completed_ids(config, "oligomer")
    if len(completed_af2_ids) < (
        len(metric_seqs) * (alphafold.LOWER_OLIGOMERS + alphafold.HIGHER_OLIGOMERS)
    ):
        af2_input_df = alphafold.make_af2_oligomer_input(metric_seqs, completed_af2_ids)
        af2_input_df.to_csv(paths.get_alphafold_input_path(config, "oligomer"))

        # run AF2 batch
        alphafold_script = alphafold.make_shell_script(config, "oligomer")
        alphafold.run_af2_batch(alphafold_script)

    # analyze ALL alphafold results
    oligomer_values_csv = Path(config["directory"]) / "oligomer_values.csv"
    if not oligomer_values_csv.is_file():
        oligomer_values = alphafold.compile_alphafold_metric_results(config)
    else:
        oligomer_values = pd.read_csv(oligomer_values_csv, index_col=0)

    metric_correlations = plotting.compute_metric_correlations(oligomer_values)
    plotting.plot_metric_correlations(
        oligomer_values, Path(config["directory"]) / "metric_plots"
    )
    print(metric_correlations)
