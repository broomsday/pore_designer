"""
Utility functions for generating sequence formats and distributions
"""

from pathlib import Path

import logomaker
import pandas as pd
import numpy as np

from designer import alphafold, proteinmpnn, sequence


def make_minimal_select_seq(
    id: str, sequence: str, oligomer: int | None = None, source: str = "wt"
) -> alphafold.SelectSeq:
    """
    Generate a SelectSeq object from a pdb_id, sequence, and oligomer count.
    """
    if oligomer is not None:
        full_sequence = [sequence] * oligomer
        unique_chains = len(sequence.split(":"))
    else:
        full_sequence = sequence
        unique_chains = len(set(full_sequence))

    merged_sequence = ":".join(full_sequence)

    return alphafold.SelectSeq(
        id=id,
        sequence=full_sequence,
        merged_sequence=merged_sequence,
        unique_chains=unique_chains,
        score=None,
        recovery=None,
        source=source,
        mutation=None,
        frequency=None,
        selection=source,
        top_plddt=None,
        mean_plddt=None,
        top_rmsd=None,
        mean_rmsd=None,
        top_hydrophobicity=None,
        top_oligomer=None,
        designed_oligomer_rank=None,
    )


def make_sequence_distribution(
    config: dict, directory: Path, pdb: Path
) -> dict[int, dict[str, float]]:
    """
    Make a sequence distribution and visualize with a sequence logo
    """
    distribution_path = directory / f"{pdb.stem}.json"

    if not distribution_path.is_file():
        sequences = [
            design.sequence[0]
            for design in proteinmpnn.get_all_sequences(config, pdb.stem)
        ]
        logo_df = logomaker.alignment_to_matrix(sequences)
        logo = logomaker.Logo(logo_df, color_scheme="weblogo_protein")
        logo.fig.savefig(directory / f"{pdb.stem}.png")
        distribution = sequence.make_consensus_frequency(sequences)
        sequence.save_distribution(distribution, distribution_path)
    else:
        distribution = sequence.load_distribution(distribution_path)

    return distribution


def make_difference_distributions(
    directory: Path, positive_distribution, negative_distribution, negative_name: str
):
    """
    Given a positive and a negative distribution,
    make a difference distribution and visualize with a sequence logo.
    """
    distribution_path = directory / f"{negative_name}.json"
    if not distribution_path.is_file():
        difference_distribution = sequence.compute_difference_distribution(
            positive_distribution, negative_distribution
        )
        difference_sequences = sequence.make_exact_sequence_list_from_distribution(
            difference_distribution
        )
        logo_df = logomaker.alignment_to_matrix(difference_sequences)
        logo = logomaker.Logo(logo_df, color_scheme="weblogo_protein")
        logo.fig.savefig(directory / f"{negative_name}.png")
        sequence.save_distribution(difference_distribution, distribution_path)
    else:
        difference_distribution = sequence.load_distribution(distribution_path)

    return difference_distribution


def score_for_negative_design(
    config: dict, negative_distributions: list, scores_df_path: Path
) -> pd.DataFrame:
    """
    Compute negative design similarity scores for all designs.
    """
    positives = proteinmpnn.get_all_sequences(config, Path(config["positive_pdb"]).stem)
    positives_df = pd.DataFrame().from_dict(
        {"sequence": [design.sequence[0] for design in positives]}
    )
    for negative_pdb, negative_distribution in zip(
        Path(config["negative_pdbs"]).glob("*.pdb"), negative_distributions
    ):
        positives_df[negative_pdb.stem] = positives_df["sequence"].apply(
            sequence.compute_blosum_similarity_by_frequency,
            args=[negative_distribution, 100],
        )

    positives_df["mean_similarity"] = positives_df.apply(
        lambda row: np.mean([value for value in row if isinstance(value, float)]),
        axis=1,
    )
    positives_df["max_similarity"] = positives_df.apply(
        lambda row: max([value for value in row if isinstance(value, float)]),
        axis=1,
    )
    positives_df.to_csv(scores_df_path)

    return positives_df


def divide_into_parts(total: int, num_parts: int) -> list[int]:
    """
    Divide an integer in roughly equal parts.
    """
    quotient, remainder = divmod(total, num_parts)
    parts = [quotient] * num_parts

    for i in range(remainder):
        parts[i] += 1

    return parts
