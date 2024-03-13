"""
Functions to assist in uring ProteinMPNN.
"""

from pathlib import Path
import itertools
from typing import NamedTuple
import subprocess
import shutil
import json

import numpy as np
import jsonlines
import biotite.structure as bts
from biotite.structure.io import load_structure
from tqdm import tqdm

from designer import paths, sequence
from designer.constants import MPNN_HIT_TYPE_ORDER


class MPNNSeq(NamedTuple):
    sequence: str
    merged_sequence: str
    unique_chains: int
    score: float
    recovery: float
    source: str
    mutation: float
    frequency: float
    selection: str


def get_chain_resnum_groups(structure: bts.AtomArray) -> dict[str, list[str]]:
    """
    Group chains by residue numbers
    """
    chain_resnum_groups = {}
    chains = bts.get_chains(structure)
    for chain in chains:
        resnums = structure[
            (structure.chain_id == chain) & (structure.atom_name == "CA")
        ].res_id
        str_resnums = "-".join([str(resnum) for resnum in resnums])
        if str_resnums in chain_resnum_groups:
            chain_resnum_groups[str_resnums].append(chain)
        else:
            chain_resnum_groups[str_resnums] = [chain]

    return chain_resnum_groups


def make_symmetry_dict(pdb_file: Path) -> dict:
    """
    Generate a ProteinMPNN symmetry dictionary for a given multimeric PDB.
    """
    structure = load_structure(pdb_file)

    # generate tied positions as ProteinMPNN likes them
    chain_resnum_groups = get_chain_resnum_groups(structure)
    tied_positions_list = []
    for resnums, chain_list in chain_resnum_groups.items():
        tied_positions_list.append(
            [
                {
                    chain: [i + 1] for chain in chain_list
                }  # residues are indexed from 1 not by res-num
                for i, _ in enumerate(resnums.split("-"))
            ]
        )
    tied_positions_list = list(itertools.chain(*tied_positions_list))
    tied_positions_dict = {pdb_file.stem: tied_positions_list}

    return tied_positions_dict


def save_symmetry_dict(symmetry_dict: dict, symmetry_file: Path) -> None:
    """
    Save the ProteinMPNN symmetry dictionary.
    """
    with open(symmetry_file, mode="wb") as fp:
        with jsonlines.Writer(fp) as writer:
            writer.write(symmetry_dict)


def rekey_proteinmpnn_dict(dict_path: Path) -> None:
    """
    Adjust the key for the PDB name in the symmetry dict to be "input",
    which we use in this code.
    """
    with jsonlines.open(dict_path) as reader:
        for object in reader:
            dict = object
            break

    for values in dict.values():
        provided_values = values
        break

    new_dict = {"input": provided_values}

    with jsonlines.open(dict_path, mode="w") as writer:
        writer.write(new_dict)


def make_shell_script(config: dict, pdb: Path | None = None) -> str:
    """
    Generate a shell script to run the remaining ProteinMPNN jobs needed.
    """
    if pdb is None:
        mpnn_folder = get_proteinmpnn_folder(config)
        pdb = config["input_pdb"]
        fixed = config["fixed_dict"]
    else:
        mpnn_folder = get_proteinmpnn_folder(config) / f"{pdb.stem}"
        fixed = None  # fixed positions not yet supported for negative design

    tied = config["symmetry_dict"]
    mpnn_folder.mkdir(exist_ok=True, parents=True)

    designs = get_num_to_design(config)
    if designs <= 0:
        raise ValueError("No designs left to do")

    shell_script = "#!/bin/bash\n\n"
    shell_script += f"source {paths.get_conda_source_path()}\n"
    shell_script += "conda deactivate\n"
    shell_script += f"conda activate {paths.get_proteinmpnn_source_path()}\n"
    shell_script += f"python {paths.get_proteinmpnn_path()}/protein_mpnn_run.py "
    shell_script += f"--pdb_path {pdb} "
    shell_script += f"--tied_positions_jsonl {tied} "
    if fixed is not None:
        shell_script += f"--fixed_positions_jsonl {fixed} "
    shell_script += f"--out_folder {mpnn_folder} "
    shell_script += f"--num_seq_per_target {designs} "
    shell_script += f"--sampling_temp {config['temperature_mpnn']} "
    shell_script += "--use_soluble_model "
    shell_script += "--omit_AAs C"

    return shell_script


def get_proteinmpnn_folder(config: dict) -> Path:
    """
    Return the top-level ProteinMPNN results folder.
    """
    return Path(config["directory"]) / "ProteinMPNN"


def run_proteinmpnn(shell_script: str) -> None:
    """
    Run ProteinMPNN.
    """
    print("Running ProteinMPNN")
    subprocess.run(
        shell_script,
        shell=True,
        executable="/bin/bash",
        stdout=subprocess.DEVNULL,
        check=False,
    )


def rename_existing_results(config: dict, pdb: str | None = None) -> None:
    """
    Rename existing results files.
    """
    mpnn_folder = get_proteinmpnn_folder(config)
    if pdb is None:
        seqs_folder = mpnn_folder / "seqs"
    else:
        seqs_folder = mpnn_folder / f"{pdb.stem}" / "seqs"

    recent_results_path = seqs_folder / "input.fa"
    if recent_results_path.is_file():
        result_count = len(list(seqs_folder.glob("*.fa")))
        new_results_path = recent_results_path.with_suffix(f".{result_count}.fa")
        shutil.move(recent_results_path, new_results_path)


def get_num_to_design(config: dict, pdb: str | None = None) -> int:
    """
    Determine how many designs still need to be made.
    """
    num_to_design = config["num_mpnn"]

    designed_seqs = get_all_sequences(config, pdb)
    num_wt_seqs = len(list((get_proteinmpnn_folder(config) / "seqs").glob("*.fa")))

    return num_to_design - (len(designed_seqs) - num_wt_seqs)


def get_sequences(fasta: Path) -> list[MPNNSeq]:
    """
    Read a ProteinMPNN output file and return a list of designs.
    """
    with open(fasta, mode="r", encoding="utf-8") as fa_file:
        lines = fa_file.readlines()

    designs = []
    for i in range(0, len(lines), 2):
        designs.append(parse_design(lines[i], lines[i + 1]))

    return designs


def get_all_sequences(config: dict, pdb: str | None = None) -> list[MPNNSeq]:
    """
    Get all designed sequences.
    """
    if pdb is None:
        mpnn_folder = get_proteinmpnn_folder(config)
    else:
        mpnn_folder = get_proteinmpnn_folder(config) / f"{pdb.stem}"

    seqs = []
    for seq_file in (mpnn_folder / "seqs").glob("*.fa"):
        seqs.extend(get_sequences(seq_file))

    return seqs


def get_wt_sequence(sequences: list[MPNNSeq]) -> MPNNSeq | None:
    """
    Pull out the WT sequence.
    """
    for sequence in sequences:
        if sequence.source == "wt":
            return sequence

    return None


def select_top_sequences(config: dict) -> list[MPNNSeq]:
    """
    Compile the top sequences designed by ProteinMPNN.
    """
    # load the sequences
    seqs = get_all_sequences(config)
    wt_sequence = get_wt_sequence(seqs)

    # get the consensus sequence and frequencies used for comparison
    frequencies = sequence.make_consensus_frequency(
        [seq.merged_sequence for seq in seqs]
    )
    merged_consensus_sequence = sequence.make_consensus_sequence(frequencies)
    consensus_sequence = merged_consensus_sequence.split(":")

    # add the consensus and frequency scores for all the sequences
    unique_seq_list = []
    scored_seqs = []
    for seq in tqdm(seqs, desc="Scoring sequences"):
        if seq.sequence in unique_seq_list:
            continue
        if seq.source == "wt":
            continue

        unique_seq_list.append(seq.sequence)
        scored_seqs.append(
            MPNNSeq(
                sequence=seq.sequence,
                merged_sequence=seq.merged_sequence,
                unique_chains=seq.unique_chains,
                score=seq.score,
                recovery=seq.recovery,
                source=seq.source,
                mutation=sequence.compute_identity(
                    merged_consensus_sequence, seq.merged_sequence
                ),
                frequency=sequence.compute_blosum_similarity_by_frequency(
                    seq.merged_sequence, frequencies
                ),
                selection=None,
            )
        )

    # start the list of selected sequences with the consensus
    selected_seqs = [
        MPNNSeq(
            sequence=consensus_sequence,
            merged_sequence=merged_consensus_sequence,
            unique_chains=len(set(consensus_sequence)),
            score=np.nan,
            recovery=sequence.compute_identity(
                wt_sequence.merged_sequence, merged_consensus_sequence
            ),
            source="consensus",
            mutation=sequence.compute_identity(
                merged_consensus_sequence, merged_consensus_sequence
            ),
            frequency=sequence.compute_blosum_similarity_by_frequency(
                merged_consensus_sequence, frequencies
            ),
            selection="consensus",
        )
    ]

    # add the top sequences by the remaining metrics
    selected_types = get_selection_types(
        config, [seq.selection for seq in selected_seqs]
    )

    num_score = selected_types.count("score")
    num_recovery = selected_types.count("recovery")
    num_mutation = selected_types.count("mutation")
    num_frequency = selected_types.count("frequency")

    scored_seqs.sort(key=lambda mpnnseq: mpnnseq.score, reverse=False)
    selected_seqs.extend(scored_seqs[:num_score])
    scored_seqs.sort(key=lambda mpnnseq: mpnnseq.recovery, reverse=True)
    selected_seqs.extend(scored_seqs[:num_recovery])
    scored_seqs.sort(key=lambda mpnnseq: mpnnseq.mutation, reverse=True)
    selected_seqs.extend(scored_seqs[:num_mutation])
    scored_seqs.sort(key=lambda mpnnseq: mpnnseq.frequency, reverse=True)
    selected_seqs.extend(scored_seqs[:num_frequency])

    return selected_seqs


def get_selection_types(config: dict, preselected: list[str]) -> list[str]:
    """
    Get a list of the different selection types (e.g. top score, top recovery)
    that will be used to build the final selection.
    """
    num_to_select = config["num_af2"] - len(preselected)
    selected_types = (
        MPNN_HIT_TYPE_ORDER * (int(num_to_select / len(MPNN_HIT_TYPE_ORDER)) + 1)
    )[:num_to_select]

    type_order = {type: i for i, type in enumerate(MPNN_HIT_TYPE_ORDER)}
    selected_types.sort(key=lambda type: type_order.get(type, np.nan))
    selected_types = preselected + selected_types

    return selected_types


def save_top_sequences(config: dict, seqs: list[MPNNSeq]) -> None:
    """
    Save the top sequences for later loading.
    """
    seqs_dict = [seq._asdict() for seq in seqs]
    with open(
        Path(config["directory"]) / "top_mpnn_results.json", mode="w", encoding="utf-8"
    ) as seq_file:
        json.dump(seqs_dict, seq_file)


def load_top_sequences(config: dict) -> list[MPNNSeq]:
    """
    Load a previously compiled set of top sequences.
    """
    with open(
        Path(config["directory"]) / "top_mpnn_results.json", mode="r", encoding="utf-8"
    ) as seq_file:
        proteinmpnn_seqs_dict = json.load(seq_file)

    # convert back to a list of MPNNSeq objects
    selection_types = get_selection_types(
        config,
        [
            seq["selection"]
            for seq in proteinmpnn_seqs_dict
            if seq["source"] != "design"
        ],
    )

    proteinmpnn_seqs = [
        MPNNSeq(
            sequence=seq["sequence"],
            merged_sequence=seq["merged_sequence"],
            unique_chains=seq["unique_chains"],
            score=seq["score"],
            recovery=seq["recovery"],
            source=seq["source"],
            mutation=seq["mutation"],
            frequency=seq["frequency"],
            selection=selection_types[i],
        )
        for i, seq in enumerate(proteinmpnn_seqs_dict)
    ]

    return proteinmpnn_seqs


def have_top_results(config: dict) -> bool:
    """
    Check if we have already parsed the top results file
    """
    folder = Path(config["directory"])
    if (folder / "top_mpnn_results.json").is_file():
        return True

    return False


def parse_design(annotation: str, sequence: str) -> MPNNSeq:
    """
    From the ProteinMPNN output lines, build the MPNNSeq data.
    """
    sequences = sequence.strip().split("/")

    unique_sequences = []
    for seq in sequences:
        if seq not in unique_sequences:
            unique_sequences.append(seq)
    unique_chains = len(set(sequences))

    merged_sequence = ":".join(sequences)

    score = [
        float(info.split("=")[1])
        for info in annotation.split(",")
        if "score" in info.split("=")[0]
    ][0]

    if "seq_recovery" in annotation:
        recovery = [
            float(info.split("=")[1])
            for info in annotation.split(",")
            if "seq_recovery" in info.split("=")[0]
        ][0]
    else:
        recovery = 1.0

    if "sample" in annotation:
        source = "design"
    else:
        source = "wt"

    return MPNNSeq(
        sequence=sequences,
        merged_sequence=merged_sequence,
        unique_chains=unique_chains,
        score=score,
        recovery=recovery,
        source=source,
        mutation=np.nan,
        frequency=np.nan,
        selection=None,
    )
