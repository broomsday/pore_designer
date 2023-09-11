"""
Functions to parse and manipulate PDBs.
"""


from pathlib import Path

import biotite.structure as bts
from biotite.structure.io import load_structure


def clean_pdb(pdb_file: Path) -> bts.AtomArray:
    """
    Load a PDB file and return a cleaned structure.
    """
    structure = load_structure(pdb_file)

    structure = structure[bts.filter_canonical_amino_acids(structure)]
    structure = structure[structure.element != "H"]

    return structure


def get_sequence(structure: bts.AtomArray) -> list[str]:
    """
    Return the triple letter sequence of a structure.
    """
    _, names = bts.get_residues(structure)
    return names


def get_multimer_state(structure: bts.AtomArray) -> int:
    """
    Return the number of identitcal chains.
    """
    total_chains = bts.get_chain_count(structure)
    chain_lengths = {len(get_sequence(chain)) for chain in bts.chain_iter(structure)}

    if len(chain_lengths) != 1:
        raise NotImplementedError("Chains have different lengths")

    return total_chains
