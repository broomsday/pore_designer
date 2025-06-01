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

from designer import paths, sequence, pdb
from designer.constants import MPNN_HIT_TYPE_ORDER, PSSM_AA_ORDER


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


def save_proteinmpnn_jsonl_dict(dict: dict, save_path: Path) -> None:
    """
    Save a ProteinMPNN jsonl dictionary.
    """
    with open(save_path, mode="wb") as fp:
        with jsonlines.Writer(fp) as writer:
            writer.write(dict)


def save_symmetry_dict(symmetry_dict: dict, symmetry_file: Path) -> None:
    """
    Save the ProteinMPNN symmetry dictionary.
    """
    save_proteinmpnn_jsonl_dict(symmetry_dict, symmetry_file)


def make_aa_bias_dict(amino_acids: list[str], bias: float) -> dict:
    """
    Make a bias dict for biasing amino acids.
    """
    return {amino_acid: bias for amino_acid in amino_acids}


def make_symmetric_fixed_positions_dict(input_stem: str, chains: list[str], positions: list[int]) -> dict:
    """
    Generate a dictionary of fixed positions for saving as a jsonl.

    Will fix ALL chains provided at the same positions.
    """
    return {
        input_stem: {
            chain: positions for chain in chains
        }
    }


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


def make_shell_script(
    config: dict,
    pdb_file: Path | None = None,
    out_folder: str | None = None,
    designs: int | None = None,
) -> str:
    """
    Generate a shell script to run the remaining ProteinMPNN jobs needed.
    """
    if pdb_file is None:
        pdb_file = config["input_pdb"]

    if out_folder is None:
        mpnn_folder = get_proteinmpnn_folder(config)
    else:
        mpnn_folder = get_proteinmpnn_folder(config) / out_folder
    mpnn_folder.mkdir(exist_ok=True, parents=True)

    if designs is None:
        designs = get_num_to_design(config)
        if designs <= 0:
            raise ValueError("No designs left to do")

    tied = config.get("symmetry_dict", None)
    fixed = config.get("fixed_dict", None)
    pssm = config.get("pssm_dict", None)
    pssm_weight = config.get("bias_weight", None)

    shell_script = "#!/bin/bash\n\n"
    shell_script += f"source {paths.get_proteinmpnn_conda_source_path()}\n"
    shell_script += "conda deactivate\n"
    shell_script += f"conda activate {paths.get_proteinmpnn_env_path()}\n"
    shell_script += f"python {paths.get_proteinmpnn_path()}/protein_mpnn_run.py "
    shell_script += f"--pdb_path {pdb_file} "
    if tied is not None:
        shell_script += f"--tied_positions_jsonl {tied} "
    if fixed is not None:
        shell_script += f"--fixed_positions_jsonl {fixed} "
    if pssm is not None:
        shell_script += f"--pssm_jsonl {pssm} "
        shell_script += "--pssm_bias_flag 1 "
        if pssm_weight is not None:
            shell_script += f"--pssm_multi {pssm_weight} "
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


def rename_existing_results(config: dict, out_folder: str | None = None) -> None:
    """
    Rename existing results files.
    """
    mpnn_folder = get_proteinmpnn_folder(config)
    if out_folder is None:
        seqs_folder = mpnn_folder / "seqs"
    else:
        seqs_folder = mpnn_folder / out_folder / "seqs"

    recent_results_path = seqs_folder / "input.fa"
    if recent_results_path.is_file():
        result_count = len(list(seqs_folder.glob("*.fa")))
        new_results_path = recent_results_path.with_suffix(f".{result_count}.fa")
        shutil.move(recent_results_path, new_results_path)


def get_num_to_design(
    config: dict, out_folder: str | None = None, num_to_design: int | None = None
) -> int:
    """
    Determine how many designs still need to be made.
    """
    if num_to_design is None:
        num_to_design = config["num_mpnn"]

    designed_seqs = get_all_sequences(config, out_folder=out_folder)
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


def get_all_sequences(config: dict, out_folder: str | None = None) -> list[MPNNSeq]:
    """
    Get all designed sequences and the WT.
    """
    if out_folder is None:
        mpnn_folder = get_proteinmpnn_folder(config)
    else:
        mpnn_folder = get_proteinmpnn_folder(config) / out_folder

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


def bias_from_counts(counts: dict[str, int], aa_order: list[str]) -> list[float]:
    """
    Take a dictionary of amino acids with their counts and return just the frequencies.

    Reorder to match the `aa_order`
    """
    total_count = sum(counts.values())
    return [counts.get(amino_acid, 0.0) / total_count for amino_acid in aa_order]


def make_pssm_from_distribution(
    config: dict, distribution: dict[int, dict[str, float]]
) -> dict[str, dict[str, dict]]:
    """
    Generate a dictionary in the format that ProteinMPNN expects for a PSSM.

    This has the following structure:
    1. A key for the PDB
    2. A key for each chain in the PDB
    3. For each chain, the following keys: pssm_coef, pssm_bias, pssm_log_odds.
        - Each of these has an entry for each residue
        - `pssm_coef` allows a custom coefficient per position, but can just be 1.0 by default
        - `pssm_bias` has a length of 21 for each residue and gives the odds per amino acid
        - `pssm_log_odds` is optional?, examples have 1.0 for each bias entry

    Example for a 2-chain, 2-residue structure:
        {
            "4JPP": {
                "A": {
                    "pssm_coef": [1.0, 1.0],
                    "pssm_bias": [
                        [0.05, 0.05, ..., 0.05. 0.0],
                        [0.05, 0.05, ..., 0.05. 0.0],
                    ],
                    "pssm_log_odds": [
                        [1.0, 1.0, ..., 1.0, 1.0],
                        [1.0, 1.0, ..., 1.0, 1.0],
                    ],
                },
                "B": {
                    "pssm_coef": [1.0, 1.0],
                    "pssm_bias": [
                        [0.05, 0.05, ..., 0.05. 0.0],
                        [0.05, 0.05, ..., 0.05. 0.0],
                    ],
                    "pssm_log_odds": [
                        [1.0, 1.0, ..., 1.0, 1.0],
                        [1.0, 1.0, ..., 1.0, 1.0],
                    ],
                },
            }
        }
    """
    # get the protein name and start the PSSM dict
    pdb_name = Path(config["positive_pdb"]).stem
    pssm = {pdb_name: {}}

    # get the chain IDs
    structure = load_structure(config["positive_pdb"])
    chains = bts.get_chains(structure)

    # generate the PSSM coef, bias, and log_odds
    pssm_coef = [1.0 for _ in distribution]
    pssm_log_odds = [[1.0 for _ in PSSM_AA_ORDER] for _ in distribution]
    pssm_bias = [
        bias_from_counts(counts, PSSM_AA_ORDER) for _, counts in distribution.items()
    ]

    # apply the distribution across each chain
    for chain in chains:
        pssm[pdb_name][chain] = {
            "pssm_coef": pssm_coef,
            "pssm_bias": pssm_bias,
            "pssm_log_odds": pssm_log_odds,
        }

    return pssm


def save_pssm(pssm_dict: dict[str, dict[str, dict]], pssm_file: Path) -> None:
    """
    Save a PSSM to a jsonlines file
    """
    with open(pssm_file, mode="wb") as fp:
        with jsonlines.Writer(fp) as writer:
            writer.write(pssm_dict)


def compute_ca_score(score_pdb: Path, out_dir: Path) -> float:
    """
    Use the Ca-only model to compute the score for a PDB
    """
    # make the Ca-only scoring shell script
    shell_script = "#!/bin/bash\n\n"
    shell_script += f"source {paths.get_proteinmpnn_conda_source_path()}\n"
    shell_script += "conda deactivate\n"
    shell_script += f"conda activate {paths.get_proteinmpnn_env_path()}\n"
    shell_script += f"python {paths.get_proteinmpnn_path()}/protein_mpnn_run.py "
    shell_script += f"--pdb_path {score_pdb} "
    shell_script += f"--out_folder {out_dir} "
    shell_script += "--ca_only "
    shell_script += "--score_only 1 "

    # check to see if the file exists already and run if not
    score_file = out_dir / "score_only" / f"{score_pdb.stem}_pdb.npz"
    if not score_file.is_file():
        # run ProteinMPNN
        print("Running ProteinMPNN")
        subprocess.run(
            shell_script,
            shell=True,
            executable="/bin/bash",
            stdout=subprocess.DEVNULL,
            check=False,
        )

    # pull out the score
    score = np.load(score_file)["score"][0]

    return float(score)


def compute_mpnn_oligomer_scores(pdb_file: Path, config: dict) -> tuple[float, float]:
    """
    Use the Ca-only model to:

    1. Compute the score for a PDB.
    2. Compute the delta between the above and when one chain is moved out 100 angstroms.
    """
    oligomer_folder = Path(config["directory"]) / "ProteinMPNN" / "oligomer"
    oligomer_folder.mkdir(exist_ok=True, parents=True)
    score = compute_ca_score(pdb_file, oligomer_folder)

    moved_folder = Path(config["directory"]) / "ProteinMPNN" / "oligomer" / "moved"
    moved_folder.mkdir(exist_ok=True)
    moved_pdb = pdb.make_chain_moved_pdb(pdb_file, moved_folder)
    moved_score = compute_ca_score(moved_pdb, moved_folder)
    delta_score = moved_score - score

    return (score, delta_score)
