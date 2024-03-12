"""
Module to hold top-level design functions
"""

from pathlib import Path
import shutil
import json

import typer
import pandas as pd

from designer import proteinmpnn, alphafold, sequence, paths, plotting, utils, pdb


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
    assert config["job_type"] == "positive"

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
        len(design_selected)
        * (config["oligomer_lower_offset"] + config["oligomer_higher_offset"])
    ):
        af2_input_df = alphafold.make_af2_oligomer_input(
            config, design_selected, completed_af2_ids
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
    assert config["job_type"] == "negative"

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

    # summarize positive and negative into distributions and save sequence logos of each
    print("Summarizing distributions and scoring sequences")
    design_summary_dir = proteinmpnn.get_proteinmpnn_folder(config) / "design"
    design_summary_dir.mkdir(exist_ok=True)

    positive_distribution = utils.make_sequence_distribution(
        config, design_summary_dir, positive_pdb
    )

    negative_distributions = []
    for negative_pdb in Path(config["negative_pdbs"]).glob("*.pdb"):
        negative_distributions.append(
            utils.make_sequence_distribution(config, design_summary_dir, negative_pdb)
        )

    # compute difference distributions for positive to each negative
    # this is positive - negative, thus, frequent amino acids maximize positive over negative
    difference_summary_dir = proteinmpnn.get_proteinmpnn_folder(config) / "difference"
    difference_summary_dir.mkdir(exist_ok=True)

    difference_distributions = []
    for negative_pdb, negative_distribution in zip(
        Path(config["negative_pdbs"]).glob("*.pdb"), negative_distributions
    ):
        difference_distributions.append(
            utils.make_difference_distributions(
                difference_summary_dir,
                positive_distribution,
                negative_distribution,
                negative_pdb.stem,
            )
        )

    # score each positive sequence by similarity to each negative distribution
    scores_df_path = Path(config["directory"]) / "scored_positives.csv"
    if not scores_df_path.is_file():
        positives_df = utils.score_for_negative_design(
            config, negative_distributions, scores_df_path
        )
    else:
        positives_df = pd.read_csv(scores_df_path)

    # choose the sequences to run AF2 on
    if (Path(config["directory"]) / "top_mpnn_results.json").is_file():
        selected_sequences = proteinmpnn.load_top_sequences(config)
    else:
        print("Selecting sequences for AF2 testing")
        num_sequences_per_type = utils.divide_into_parts(int(config["num_af2"]), 3)

        # choose 1/3 of designs based on overall similarity
        positives_df = positives_df.sort_values(by="mean_similarity")
        selected_sequences = [
            utils.make_minimal_select_seq(
                i, positives_df.iloc[i].sequence, config["multimer"], "mean_similarity"
            )
            for i in range(num_sequences_per_type[0])
        ]

        # choose 1/3 of designs based on max similarity
        positives_df = positives_df.sort_values(by="max_similarity")
        selected_sequences.extend(
            [
                utils.make_minimal_select_seq(
                    i,
                    positives_df.iloc[i].sequence,
                    config["multimer"],
                    "max_similarity",
                )
                for i in range(num_sequences_per_type[1])
            ]
        )

        # randomly sample an equal number of sequences from the difference distribution of the +1 oligomer
        for negative_pdb in Path(config["negative_pdbs"]).glob("*.pdb"):
            multimer_state = pdb.get_multimer_state(pdb.load_pdb(negative_pdb))
            if multimer_state == int(config["multimer"]) + 1:
                difference_distribution = sequence.load_distribution(
                    proteinmpnn.get_proteinmpnn_folder(config)
                    / "difference"
                    / f"{negative_pdb.stem}.json"
                )

        difference_sequences = sequence.sample_sequences_from_distribution(
            difference_distribution, num_sequences_per_type[2]
        )
        selected_sequences.extend(
            [
                utils.make_minimal_select_seq(
                    i, difference_sequences[i], config["multimer"], "difference"
                )
                for i in range(num_sequences_per_type[1])
            ]
        )

        # save the selected sequences as though they had been made by ProteinMPNN so we can use them downstream
        proteinmpnn.save_top_sequences(config, selected_sequences)

    # setup AF2 input file
    completed_af2_ids = alphafold.get_completed_ids(config, "design")
    if len(completed_af2_ids) < config["num_af2"]:
        af2_input_df = alphafold.make_af2_design_input(
            selected_sequences, completed_af2_ids
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
        len(design_selected)
        * (config["oligomer_lower_offset"] + config["oligomer_higher_offset"])
    ):
        af2_input_df = alphafold.make_af2_oligomer_input(
            config, design_selected, completed_af2_ids
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


def design_pore_mutations(
    config_path: Path = typer.Argument(
        ..., help="Pore Design config file for mutation design"
    )
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job for mutation design.
    """
    with open(config_path, mode="r", encoding="utf-8") as config_file:
        config = json.load(config_file)
    assert config["job_type"] == "mutation"

    # get the WT sequence and check that it works for our use case
    wt_structure = pdb.load_pdb(Path(config["input_pdb"]))
    wt_sequence = pdb.get_by_chain_sequence(wt_structure)
    sequence_length = list(set([len(seq) for seq in wt_sequence]))

    if len(sequence_length) > 1:
        raise NotImplementedError(
            "Mutation design not implemented for hetero-oligomers"
        )

    # if we're using a mutation distribution load it now otherwise generate a uniform distribution
    if config["mutation_distribution"] is not None:
        mutation_distribution = sequence.load_distribution(
            Path(config["mutation_distribution"])
        )
        if len(mutation_distribution) != sequence_length[0]:
            assert ValueError(
                "Provided mutation distribution's length does not match the input"
            )
    else:
        mutation_distribution = sequence.mutation_distribution_from_length(
            sequence_length[0], bias=None
        )

    # build the list of new sequences by choosing mutations
    selected_sequences = [
        utils.make_minimal_select_seq(
            i,
            sequence.mutate_sequence_by_distribution(
                wt_sequence[0], mutation_distribution, config["mutations_per_seq"]
            ),
            config["multimer"],
            "mutation",
        )
        for i in range(config["num_af2"])
    ]
    # save the mutated sequences as though they had been made by ProteinMPNN so we can use them downstream
    proteinmpnn.save_top_sequences(config, selected_sequences)

    # setup AF2 input file
    completed_af2_ids = alphafold.get_completed_ids(config, "design")
    if len(completed_af2_ids) < config["num_af2"]:
        af2_input_df = alphafold.make_af2_design_input(
            selected_sequences, completed_af2_ids
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
        len(design_selected)
        * (config["oligomer_lower_offset"] + config["oligomer_higher_offset"])
    ):
        af2_input_df = alphafold.make_af2_oligomer_input(
            config, design_selected, completed_af2_ids
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
    assert config["job_type"] == "metric"

    # build the input data into SelectSeq objects
    input_data = pd.read_csv(config["input_csv"], index_col="pdb")
    metric_seqs = [
        utils.make_minimal_select_seq(
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
        len(metric_seqs)
        * (config["oligomer_lower_offset"] + config["oligomer_higher_offset"])
    ):
        af2_input_df = alphafold.make_af2_oligomer_input(
            config, metric_seqs, completed_af2_ids
        )
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
