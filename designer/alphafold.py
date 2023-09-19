"""
Functions for using AlphaFold2.
"""


from pathlib import Path
import subprocess
from typing import NamedTuple

import pandas as pd

from designer.proteinmpnn import MPNNSeq
from designer import paths


class AFSeq(NamedTuple):
    sequence: str
    merged_sequence: str
    unique_chains: int
    score: float
    recovery: float
    source: str
    consensus: float
    frequency: float
    selection: str


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


def make_af2_oligomer_input(selected_seqs: list) -> pd.DataFrame:
    """
    Save sequences in oligomer format for checking.
    """
    # TODO: code
    pass


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


def select_top(
    config: dict,
    plddt: float = 0.9,
    rmsd: float = 1.0,
    oligomer: int | None = None,
    identity: float = 0.9,
) -> list[AFSeq]:
    """
    Select the top AF2 sequences based on a set of criteria.
    """
    # TODO: load the ProteinMPNN input seqs so that we have MPNNSeq data (proteinmpnn.load_top_sequences)

    # TODO: unzip (utils.)
    # TODO: pull out plddt (_.)
    # TODO: compute RMSD (pdb.)
    # TODO: optionally compute the oligomer check (_.)
    # TODO:
    pass


def have_selected(config: dict, phase: str) -> bool:
    """
    Check if we have analyzed and selected the best sequences.
    """
    # TODO: code
    return False


def save_selected(config: dict, phase: str, selected: list[AFSeq]) -> None:
    """
    Save the selected sequences.
    """
    # TODO: code
    pass


def load_selected(config: dict, phase: str) -> list[AFSeq]:
    """
    Load the selected sequences.
    """
    # TODO: code
    pass


def report_selected(config: dict, selected: list) -> None:
    """
    Save a dataframe of the selected sequences and copy over the AF2 PDBs.
    Also produce and copy over the oligomer-check graphs if these are multimers.
    """
    # TODO: code
    pass
