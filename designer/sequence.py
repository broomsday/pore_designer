"""
Functions for working with and computing differences between sequences.
"""

import numpy as np

from designer.constants import BLOSUM_MATRIX, CANONICAL_AMINO_ACIDS


def compute_blosum_similarity(seq_one: str, seq_two: str) -> float:
    """
    Compute the similarity between two sequences by using the blosum matrix.

    Sequences must be identical length and taken as being aligned.
    """
    assert len(seq_one) == len(seq_two)

    return np.mean(
        [BLOSUM_MATRIX[aa_one][aa_two] for aa_one, aa_two in zip(seq_one, seq_two)]
    )


def compute_blosum_similarity_by_frequency(
    sequence: str, frequency: dict[int, dict[str, float]]
) -> float:
    """
    Return the BLOSUM similarity weighted by the amino acid frequencies of the consensus.

    NOTE: in theory this should only be ~20x slower than `compute_blosum_similarity` but it is even slower...
    """
    similarity = 0.0
    for i, aa in enumerate(sequence):
        similarity += np.sum(
            [
                compute_blosum_similarity(aa, freq_aa) * frequency[i][freq_aa]
                for freq_aa in frequency[i]
            ]
        )

    return similarity


def compute_amino_acid_frequencies(amino_acids: list[str]) -> dict[str, int]:
    """
    Given a list of amino acids, return a dictionary counting each unique amino acid.
    """
    return {
        amino_acid: amino_acids.count(amino_acid) / len(amino_acids)
        for amino_acid in CANONICAL_AMINO_ACIDS
    }


def make_consensus_frequency(sequences: list[str]) -> dict[int, dict[str, float]]:
    """
    Given an aligned set of sequences, determine frequency of amino acids at each position.
    """
    return {
        position: compute_amino_acid_frequencies(
            [sequence[position] for sequence in sequences]
        )
        for position in range(len(sequences[0]))
    }


def get_highest_frequency_amino_acid(frequencies: dict[str, float]) -> str:
    """
    Given frequencies of some amino acids, return the highest one.
    """
    highest_frequency = max(frequencies.values())
    return {frequency: amino_acid for amino_acid, frequency in frequencies.items()}[
        highest_frequency
    ]


def make_consensus_sequence(consensus_frequency: dict[int, dict[str, float]]) -> str:
    """
    Given amino acid frequencies at each position, return the consensus as the higest frequency amino
    acid at each position.
    """
    return "".join(
        [
            get_highest_frequency_amino_acid(frequencies)
            for frequencies in consensus_frequency.values()
        ]
    )


def compute_identity(reference: str, query: str) -> float:
    """
    Return the fraction of identical amino acids assuming same length sequences.
    """
    assert len(reference) == len(query)
    return np.mean(
        [1 if reference[i] == query_aa else 0 for i, query_aa in enumerate(query)]
    )
