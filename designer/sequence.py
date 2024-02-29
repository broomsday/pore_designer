"""
Functions for working with and computing differences between sequences.
"""

from pathlib import Path
import json
import random

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
    sequence: str, frequency: dict[int, dict[str, float]], frequency_scaler: float = 1
) -> float:
    """
    Return the BLOSUM similarity to a weighted frequency (e.g. a consensus)
    """
    similarity = 0.0
    for i, aa in enumerate(sequence):
        similarity += np.sum(
            [
                compute_blosum_similarity(aa, freq_aa)
                * (frequency[i][freq_aa] / frequency_scaler)
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


def save_distribution(
    distribution: dict[int, dict[str, float]], out_path: Path
) -> None:
    """
    Save a sequence distribution (logo or matrix) as a json file.
    """
    with open(out_path, mode="w", encoding="utf8") as out_file:
        json.dump(distribution, out_file)


def load_distribution(in_path: Path) -> dict[int, dict[str, float]]:
    """
    Load a sequence distribution (logo or matrix) from a json file.
    """
    with open(in_path, mode="r", encoding="utf8") as in_file:
        distribution = json.load(in_file)

    # JSON serializes with keys as strings not ints, so convert back
    distribution = {
        int(position): frequencies for position, frequencies in distribution.items()
    }

    return distribution


def compute_difference_distribution(
    positive_distribution: dict[int, dict[str, float]],
    negative_distribution: dict[int, dict[str, float]],
    zero_floor: bool = True,
) -> dict[int, dict[str, float]]:
    """
    Given two distributions compute the difference.

    If `zero_floor`, negative frequencies are set to a floor of 0.
    """
    assert positive_distribution.keys() == negative_distribution.keys()

    difference_distribution = {}
    for position in positive_distribution:
        # if the two positions are equal don't bother computing a difference
        if positive_distribution[position] == negative_distribution[position]:
            difference_distribution[position] = positive_distribution[position]
            continue

        # compute the difference
        if zero_floor:
            difference = {
                amino_acid: (
                    max(
                        positive_distribution[position][amino_acid]
                        - negative_distribution[position][amino_acid],
                        0,
                    )
                )
                for amino_acid in positive_distribution[position]
            }
        else:
            difference = {
                amino_acid: (
                    positive_distribution[position][amino_acid]
                    - negative_distribution[position][amino_acid]
                )
                for amino_acid in positive_distribution[position]
            }

        # normalize after computation
        frequency_sum = sum(difference.values())
        difference = {
            amino_acid: frequency / frequency_sum
            for amino_acid, frequency in difference.items()
        }

        # add to the final distribution
        difference_distribution[position] = difference

    return difference_distribution


def normalize_to_integer_sum(
    values: list[float | int], integer_sum: int = 100
) -> list[int]:
    """
    Take a list of non-negative numbers and normalize to have a given integer sum, such that
    all values are themselves non-negative integers
    """
    # Step 1: Normalize to sum of 100 in float
    total_sum = np.sum(values)
    normalized_floats = [(x / total_sum) * integer_sum for x in values]

    # Step 2: Round to nearest integer
    rounded_ints = [round(x) for x in normalized_floats]

    # Step 3: Adjust the sum to exactly 100
    adjustment = integer_sum - sum(rounded_ints)

    # Sort items by the decimal part of their original normalized values, which were lost during rounding
    # This helps in deciding which numbers to adjust
    decimals = [x - int(x) for x in normalized_floats]
    sorted_indices = sorted(
        range(len(values)),
        key=lambda i: decimals[i],
        reverse=True if adjustment > 0 else False,
    )

    # Adjust the rounded integers
    adjustment_count = 0
    loop_count = 0
    while adjustment_count < abs(adjustment):
        loop_count += 1
        # If adjustment is positive, increment; if negative, decrement unless that would make it negative
        index_to_adjust = sorted_indices[loop_count % len(values)]
        if adjustment > 0:
            rounded_ints[index_to_adjust] += 1
            adjustment_count += 1
        else:
            if rounded_ints[index_to_adjust] != 0:
                rounded_ints[index_to_adjust] -= 1
                adjustment_count += 1

        if loop_count > 100:
            raise ValueError

    return rounded_ints


def make_exact_sequence_list_from_distribution(
    distribution: dict[int, dict[str, float]]
) -> list[str]:
    """
    Round the amino acid distribution at each position to 1% intervals.

    Then produce a list of 100 sequences that perfectly fit the above distribution.
    """
    # normalize the distribution to a sum of 100 at each position with integer frequencies only
    for position in distribution:
        new_frequencies = normalize_to_integer_sum(
            list(distribution[position].values())
        )
        distribution[position] = {
            amino_acid: new_frequencies[i]
            for i, amino_acid in enumerate(distribution[position])
        }

    # build the list of sequences one position at a time
    sequences = ["" for _ in range(100)]
    for position in distribution:
        assert np.sum(list(distribution[position].values())) == 100

        position_sum = 0
        added_positions = 0
        for amino_acid, count in distribution[position].items():
            if count == 0:
                continue
            for i in range(count):
                sequences[i + position_sum] += amino_acid
                added_positions += 1
            position_sum += count

    return sequences


def sample_amino_acid_by_weight(amino_acid_weights: dict[str, int | float]) -> str:
    """
    Given a dictionary of amino acids as keys and weights as values,
    return a weighted random amino acid.
    """
    return random.choices(
        list(amino_acid_weights.keys()), list(amino_acid_weights.values())
    )[0]


def sample_sequences_from_distribution(
    distribution: dict[int, dict[str, float]], num_samples: int
) -> list[str]:
    """
    Given a distribution, randomly sample sequences from it.
    """
    sequences = [
        "".join(
            [
                sample_amino_acid_by_weight(amino_acid_weights)
                for _, amino_acid_weights in distribution.items()
            ]
        )
        for _ in range(num_samples)
    ]

    return sequences
