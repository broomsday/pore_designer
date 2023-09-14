"""
Functions to assist in uring ProteinMPNN.
"""


from pathlib import Path
import itertools

import jsonlines
import biotite.structure as bts
from biotite.structure.io import load_structure


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
                # {chain.id: [int(resnum)] for chain in chain_list}
                {
                    chain: [i + 1] for chain in chain_list
                }  # residues are indexed from 1 not by res-num
                for i, resnum in enumerate(resnums.split("-"))
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
