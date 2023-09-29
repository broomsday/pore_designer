"""
Functions for using AlphaFold2.
"""


from pathlib import Path
import subprocess
import shutil
from typing import NamedTuple
import json

import pandas as pd
from tqdm import tqdm

from designer.proteinmpnn import MPNNSeq
from designer import proteinmpnn, file_utils, pdb, paths, sequence


LOWER_OLIGOMERS = 3
HIGHER_OLIGOMERS = 4


class SelectSeq(NamedTuple):
    id: str
    sequence: str
    merged_sequence: str
    unique_chains: int
    score: float
    recovery: float
    source: str
    mutation: float
    frequency: float
    selection: str
    top_plddt: float
    mean_plddt: float
    top_rmsd: float
    mean_rmsd: float
    top_oligomer: int
    designed_oligomer_rank: int


def make_id(seq: MPNNSeq, index: int) -> str:
    """
    Given sequence annotations, construct an ID.
    """
    return f"{seq.selection}_{index}"


def get_completed_ids(config: dict, phase: str) -> list[str]:
    """
    Get a list of all sequences that have already been predicted by AF2.
    """
    alphafold_dir = Path(config["directory"]) / "AlphaFold" / phase / "outputs"

    return [
        result.name.replace(".result.zip", "")
        for result in alphafold_dir.glob("*.result.zip")
    ]


def make_af2_design_input(
    seqs: list[MPNNSeq], completed_ids: list[str]
) -> pd.DataFrame:
    """
    Save the list of sequences as AF2 input.
    """
    seqs_dict = {
        "id": [make_id(seq, i) for i, seq in enumerate(seqs)],
        "sequence": [":".join(seq.sequence) for seq in seqs],
    }
    af2_inputs = pd.DataFrame().from_dict(seqs_dict)
    af2_inputs = af2_inputs[~af2_inputs["id"].isin(completed_ids)]

    af2_inputs = af2_inputs.set_index("id")
    af2_inputs.index.name = "id"

    return af2_inputs


def make_af2_oligomer_input(
    selected_seqs: list, completed_ids: list[str]
) -> pd.DataFrame:
    """
    Save sequences in oligomer format for checking.
    """
    oligomer_dict = {"id": [], "sequence": []}
    for selected_seq in selected_seqs:
        designed_oligomer = len(selected_seq.sequence)
        subunit_seq = [
            selected_seq.sequence[i] for i in range(selected_seq.unique_chains)
        ]

        for oligomer_offset in range(LOWER_OLIGOMERS):
            oligomers = designed_oligomer - oligomer_offset
            oligomer_dict["id"].append(f"{selected_seq.id}_{oligomers}")
            oligomer_dict["sequence"].append(":".join(subunit_seq * oligomers))

        for oligomer_offset in range(HIGHER_OLIGOMERS):
            oligomers = designed_oligomer + oligomer_offset
            oligomer_dict["id"].append(f"{selected_seq.id}_{oligomers}")
            oligomer_dict["sequence"].append(":".join(subunit_seq * oligomers))

    af2_inputs = pd.DataFrame().from_dict(oligomer_dict)
    af2_inputs = af2_inputs[~af2_inputs["id"].isin(completed_ids)]

    af2_inputs = af2_inputs.set_index("id")
    af2_inputs.index.name = "id"

    return af2_inputs


def make_shell_script(config: dict, phase: str) -> str:
    """
    Make the shell script for AF2 running.
    """
    af2_run_path = Path(config["directory"]) / "AlphaFold" / f"{phase}" / "outputs"
    af2_run_path.mkdir(exist_ok=True, parents=True)

    if config["multimer"] == 1:
        model_type = "alphafold2_ptm"
    else:
        model_type = "alphafold2_multimer_v3"

    shell_script = "#!/bin/bash\n\n"
    shell_script = f"cd {af2_run_path.absolute()}\n"
    shell_script += f"source {paths.get_conda_source_path()}\n"
    shell_script += "conda deactivate\n"
    shell_script += f"conda activate {paths.get_alphafold_source_path()}\n"
    shell_script += f"colabfold_batch {paths.get_alphafold_input_path(config, phase)} {af2_run_path} "
    shell_script += "--recompile-padding 1 "
    shell_script += "--msa-mode single_sequence "
    shell_script += f"--model-type {model_type} "
    shell_script += f"--num-recycle {config['recycle_af2']} "
    shell_script += "--num-models 5 "
    shell_script += "--zip"

    return shell_script


def run_af2_batch(shell_script: str) -> None:
    """
    Execute AF2 script.
    """
    print("Running AlphaFold2")
    subprocess.run(
        shell_script,
        shell=True,
        executable="/bin/bash",
        # stdout=subprocess.DEVNULL,
        check=False,
    )


def compile_alphafold_results(config: dict, phase: str) -> list[SelectSeq]:
    """
    Compile all the alphafold and proteinmpnn metrics together.
    """
    # load the ProteinMPNN sequence data
    proteinmpnn_seqs = proteinmpnn.load_top_sequences(config)
    # build a list of just the merged sequences for cross-reference lookup
    proteinmpnn_merged_sequences = [seq.merged_sequence for seq in proteinmpnn_seqs]
    # load the FULL alphafold input for cross-referencing with above
    alphafold_input = make_af2_design_input(proteinmpnn_seqs, [])

    # cleanup the alphafold results (last run with have left unzipped contents)
    alphafold_result_dir = Path(config["directory"]) / "AlphaFold" / phase / "outputs"
    file_utils.keep_files_by_suffix(alphafold_result_dir, [".zip"])

    # for each result, pull out alphafold metrics and combine with proteinmpnn metrics
    seqs = []
    for result_file in tqdm(
        list(alphafold_result_dir.glob("*.zip")), desc="Computing Alphafold Results"
    ):
        sequence_id = result_file.stem.replace(".result", "")
        result_dir = alphafold_result_dir / sequence_id

        # unzip result file
        file_utils.unzip_af2_prediction(
            result_file,
            result_dir,
        )

        # pull out plddts
        top_plddt = file_utils.get_average_plddt(result_dir, top_only=True)
        mean_plddt = file_utils.get_average_plddt(result_dir, top_only=False)

        # compute RMSD
        input_pdb = paths.get_input_pdb_path(config)
        top_rmsd = pdb.compute_rmsd_to_template(result_dir, input_pdb, top_only=True)
        mean_rmsd = pdb.compute_rmsd_to_template(result_dir, input_pdb, top_only=False)

        # optionally compute the oligomer check
        if (config["multimer"] > 1) and (phase == "oligomer"):
            oligomer_ranks = compute_oligomer_ranks(config)
            top_oligomer = min(oligomer_ranks, key=oligomer_ranks.get)
            designed_oligomer_rank = oligomer_ranks[config["multimer"]]
        else:
            top_oligomer = None
            designed_oligomer_rank = None

        # cross-reference the sequence with proteinmpnn to get proteinmpnn metrics
        merged_sequence = alphafold_input.loc[sequence_id].sequence.replace(":", "_")
        sequence_index = proteinmpnn_merged_sequences.index(merged_sequence)
        proteinmpnn_seq = proteinmpnn_seqs[sequence_index]

        # contruct the SelectSeq object and append
        seq = SelectSeq(
            id=sequence_id,
            sequence=proteinmpnn_seq.sequence,
            merged_sequence=proteinmpnn_seq.merged_sequence,
            unique_chains=proteinmpnn_seq.unique_chains,
            score=proteinmpnn_seq.score,
            recovery=proteinmpnn_seq.recovery,
            source=proteinmpnn_seq.source,
            mutation=proteinmpnn_seq.mutation,
            frequency=proteinmpnn_seq.frequency,
            selection=proteinmpnn_seq.selection,
            top_plddt=top_plddt,
            mean_plddt=mean_plddt,
            top_rmsd=top_rmsd,
            mean_rmsd=mean_rmsd,
            top_oligomer=top_oligomer,
            designed_oligomer_rank=designed_oligomer_rank,
        )
        seqs.append(seq)

    return seqs


def passes_criteria(
    seq: SelectSeq,
    top_plddt: float,
    mean_plddt: float,
    top_rmsd: float,
    mean_rmsd: float,
    oligomer: int,
) -> bool:
    """
    Check if the given sequence passes the AlphaFold criteria
    """
    if seq.top_plddt < top_plddt:
        return False
    if seq.mean_plddt < mean_plddt:
        return False
    if seq.top_rmsd > top_rmsd:
        return False
    if seq.mean_rmsd > mean_rmsd:
        return False
    if (seq.designed_oligomer_rank is not None) and (
        seq.designed_oligomer_rank < oligomer
    ):
        return False

    return True


def select_top(
    evaluated_seqs: list[SelectSeq],
    oligomer_results: list[SelectSeq] | None = None,
    top_plddt: float = 90,
    mean_plddt: float = 80,
    top_rmsd: float = 1.0,
    mean_rmsd: float = 2.0,
    oligomer: int | None = None,
    max_identity: float = 0.9,
) -> list[SelectSeq]:
    """
    Select the top AF2 sequences based on a set of criteria.
    """
    if oligomer_results is not None:
        raise NotImplementedError("implement oligomer check stuff")

    # filter based on plddt and RMSD
    filtered_seqs = [
        seq
        for seq in evaluated_seqs
        if passes_criteria(seq, top_plddt, mean_plddt, top_rmsd, mean_rmsd, oligomer)
    ]
    if len(filtered_seqs) == 0:
        raise ValueError(
            "plddt and/or rmsd cutoffs are too strict, no sequences found.\n"
            f"max top plddt: {max([seq.top_plddt for seq in evaluated_seqs])}\n"
            f"max mean plddt: {max([seq.mean_plddt for seq in evaluated_seqs])}\n"
            f"min top rmsd: {min([seq.top_rmsd for seq in evaluated_seqs])}\n"
            f"min mean rmsd: {min([seq.mean_rmsd for seq in evaluated_seqs])}\n"
        )

    # keep only the best given the identity cutoff
    filtered_seqs.sort(key=lambda seq: seq.top_rmsd)
    selected_seqs = [filtered_seqs.pop(0)]
    for filtered_seq in filtered_seqs:
        fails_identity = False
        for selected_seq in selected_seqs:
            identity = sequence.compute_identity(
                filtered_seq.merged_sequence, selected_seq.merged_sequence
            )
            if identity > max_identity:
                fails_identity = True
                break
        if not fails_identity:
            selected_seqs.append(filtered_seq)

    return selected_seqs


def have_alphafold(config: dict, phase: str, stage: str) -> bool:
    """
    Check if we have analyzed the Alphafold runs.
    """
    if stage == "results":
        alphafold_file = paths.get_alphafold_results_path(config, phase)
    elif stage == "selected":
        alphafold_file = paths.get_alphafold_selected_path(config, phase)

    return alphafold_file.is_file()


def save_alphafold(config: dict, phase: str, stage: str, seqs: list[SelectSeq]) -> None:
    """
    Save the selected sequences.
    """
    if stage == "results":
        alphafold_file = paths.get_alphafold_results_path(config, phase)
    elif stage == "selected":
        alphafold_file = paths.get_alphafold_selected_path(config, phase)

    seqs_dict = [seq._asdict() for seq in seqs]
    with open(
        alphafold_file,
        mode="w",
        encoding="utf-8",
    ) as seq_file:
        json.dump(seqs_dict, seq_file)


def load_alphafold(config: dict, phase: str, stage: str) -> list[SelectSeq]:
    """
    Load the selected sequences.
    """
    if stage == "results":
        alphafold_file = paths.get_alphafold_results_path(config, phase)
    elif stage == "selected":
        alphafold_file = paths.get_alphafold_selected_path(config, phase)

    with open(
        alphafold_file,
        mode="r",
        encoding="utf-8",
    ) as seq_file:
        seq_dict = json.load(seq_file)

    selected_seqs = [
        SelectSeq(
            id=seq["id"],
            sequence=seq["sequence"],
            merged_sequence=seq["merged_sequence"],
            unique_chains=seq["unique_chains"],
            score=seq["score"],
            recovery=seq["recovery"],
            source=seq["source"],
            mutation=seq["mutation"],
            frequency=seq["frequency"],
            selection=seq["selection"],
            top_plddt=seq["top_plddt"],
            mean_plddt=seq["mean_plddt"],
            top_rmsd=seq["top_rmsd"],
            mean_rmsd=seq["mean_rmsd"],
            top_oligomer=seq["top_oligomer"],
            designed_oligomer_rank=seq["designed_oligomer_rank"],
        )
        for seq in seq_dict
    ]

    return selected_seqs


def compute_oligomer_ranks(config: dict) -> dict[int:float]:
    """
    Analyze the oligomer check data and report the plddt for each oligomer.
    """
    pass
    raise NotImplementedError("code oligomer rank stuff once we have data")


def report_selected(config: dict, selected: list[SelectSeq]) -> None:
    """
    Print out selected sequence report and move relevant files into a final report dir.
    Save the report as a dataframe also.
    """
    # make the final directory
    final_directory = Path(config["directory"]) / "final_selected"
    final_directory.mkdir(exist_ok=True)

    for seq in selected:
        seq_directory = final_directory / str(seq.id)
        seq_directory.mkdir(exist_ok=True)

        # copy over the PDBs for each winner
        results_dir = (
            Path(config["directory"]) / "AlphaFold" / "design" / "outputs" / seq.id
        )
        top_pdb = list(results_dir.glob("*rank_001*.pdb"))[0]
        shutil.copy(top_pdb, seq_directory / f"{seq.id}.pdb")

    # if multimer, make the oligomer plots
    if config["multimer"] > 1:
        plot_oligomer_check(config)
        pass

    # build a dataframe of the results, save, and print
    selected_dict = [seq._asdict() for seq in selected]
    selected_df = pd.DataFrame().from_dict(selected_dict)
    selected_df.to_csv(Path(config["directory"]) / "final_selected.csv")
    print(selected_df)


def plot_oligomer_check(config) -> None:
    """
    Load the plddt values for each oligomer and plot
    """
    raise NotImplementedError("code oligomer plotting")
