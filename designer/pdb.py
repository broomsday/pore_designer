"""
Functions to parse and manipulate PDBs.
"""

from pathlib import Path
import itertools

import numpy as np
import biotite.structure as bts
from biotite.structure.io import load_structure, save_structure

from designer.constants import AMINO_ACID_THREE_TO_ONE, HYDROPHOBIC_ELEMENTS


def load_pdb(pdb_file: Path) -> bts.AtomArray:
    """
    Load a PDB as a biotite AtomArray
    """
    return load_structure(pdb_file)


def clean_pdb(pdb_file: Path) -> bts.AtomArray:
    """
    Load a PDB file and return a cleaned structure.
    """
    structure = load_structure(pdb_file)

    structure = structure[bts.filter_canonical_amino_acids(structure)]
    structure = structure[structure.element != "H"]

    return structure


def get_sequence(structure: bts.AtomArray, mode: str = "triple") -> list[str]:
    """
    Return the single or triple (default) letter sequence of a structure.
    """
    _, names = bts.get_residues(structure)

    if mode == "triple":
        return names

    names = [AMINO_ACID_THREE_TO_ONE[name] for name in names]
    return names


def get_by_chain_sequence(structure: bts.AtomArray) -> list[list[str]]:
    """
    Return a list of the single letter sequence for each chain
    """
    return [
        "".join(get_sequence(chain, mode="single"))
        for chain in bts.chain_iter(structure)
    ]


def get_multimer_state(structure: bts.AtomArray) -> int:
    """
    Return the number of identitcal chains.
    """
    total_chains = bts.get_chain_count(structure)
    chain_lengths = {len(get_sequence(chain)) for chain in bts.chain_iter(structure)}

    if len(chain_lengths) != 1:
        raise NotImplementedError("Chains have different lengths")

    return total_chains


def reorder_structure_by_chains(
    structure: bts.AtomArray, chain_order: list[str]
) -> bts.AtomArray:
    """
    Given a structure and a list of chain ids, reorder the atoms in the structure to match the order of chains
    in the input list.
    """
    # start with the first chain
    reordered_structure = structure[np.isin(structure.chain_id, chain_order[0])]

    # if there are additional chains, add them on in order
    if len(chain_order) > 1:
        for chain in chain_order[1:]:
            reordered_structure += structure[np.isin(structure.chain_id, chain)]

    return reordered_structure


def get_next_best_chains(
    template: bts.AtomArray,
    predicted: bts.AtomArray,
    ordered_template_chains: list[str],
    ordered_predicted_chains: list[str],
    unordered_template_chains: set[str],
    unordered_predicted_chains: set[str],
) -> tuple[str, str]:
    """
    Given a template and predicted structure and the ordered list of current matching chains, find the next best
    matching chains that minimize the RMSD between structures.
    """
    # ordered the structures based on our current progress (if there is any)
    if (len(ordered_template_chains) >= 1) and (len(ordered_predicted_chains) >= 1):
        reordered_template = reorder_structure_by_chains(
            template, ordered_template_chains
        )
        reordered_predicted = reorder_structure_by_chains(
            predicted, ordered_predicted_chains
        )
    else:
        reordered_template = None
        reordered_predicted = None

    # now find the next best chain
    best_case = {"template_chain": "", "predicted_chain": "", "rmsd": np.inf}
    for template_chain, predicted_chain in itertools.product(
        unordered_template_chains, unordered_predicted_chains
    ):
        # temporarily add these chains
        if (reordered_template is not None) and (reordered_predicted is not None):
            test_template = (
                reordered_template
                + template[np.isin(template.chain_id, template_chain)]
            )
            test_predicted = (
                reordered_predicted
                + predicted[np.isin(predicted.chain_id, predicted_chain)]
            )
        else:
            test_template = template[np.isin(template.chain_id, template_chain)]
            test_predicted = predicted[np.isin(predicted.chain_id, predicted_chain)]

        # compute the RMSD and check vs. our best RMSD
        test_superimposed, _ = bts.superimpose(test_template, test_predicted)
        rmsd = bts.rmsd(test_template, test_superimposed)
        if rmsd < best_case["rmsd"]:
            best_case["rmsd"] = rmsd
            best_case["template_chain"] = template_chain
            best_case["predicted_chain"] = predicted_chain

    return best_case["template_chain"], best_case["predicted_chain"]


def iteratively_superimpose_chains(
    template: bts.AtomArray, predicted: bts.AtomArray
) -> tuple[bts.AtomArray, bts.AtomArray]:
    """
    Superimpose and compute RMSD for each chain permutation between template and predicted.
    Identify the lowest RMSD chain pair and superimpose on those.
    Repeat with the remaining chains.
    When no chains remain, return the superimposed structure.
    """
    unordered_template_chains = set(template.chain_id)
    unordered_predicted_chains = set(predicted.chain_id)
    ordered_template_chains = []
    ordered_predicted_chains = []

    # progressively order/matche the chains to minimze the RMSD until all chains have been matched
    while (len(unordered_template_chains) > 0) or (len(unordered_predicted_chains) > 0):
        template_chain, predicted_chain = get_next_best_chains(
            template,
            predicted,
            ordered_template_chains,
            ordered_predicted_chains,
            unordered_template_chains,
            unordered_predicted_chains,
        )
        unordered_template_chains.remove(template_chain)
        unordered_predicted_chains.remove(predicted_chain)
        ordered_template_chains.append(template_chain)
        ordered_predicted_chains.append(predicted_chain)

    # reorder the structure now that we have the correct chain order
    ordered_template = reorder_structure_by_chains(template, ordered_template_chains)
    ordered_predicted = reorder_structure_by_chains(predicted, ordered_predicted_chains)

    # align
    superimposed_ordered_predicted, _ = bts.superimpose(
        ordered_template, ordered_predicted
    )
    return ordered_template, superimposed_ordered_predicted


def compute_rmsd_to_template(
    predicted_folder: Path, template_pdb: Path, top_only: bool = True
) -> float:
    """
    Compute the Ca RMSD between the AF2 prediction and the reference.

    Note this calculation is somewhat expensive because the order of chains in the template may not match the order
    in the AF2 prediction, thus a trivial 1:1 correspondance of atoms cannot be assumed and we need to reorder the
    chains into the lowest RMSD configuration.
    """
    # load the template and reduce to Ca only
    template = load_structure(template_pdb)
    template = template[np.isin(template.atom_name, "CA")]

    # load the predicted and reduce to Ca only
    predicted_pdbs = list(predicted_folder.glob("*.pdb"))
    if top_only:
        predicted_pdbs = [
            predicted_pdb
            for predicted_pdb in predicted_pdbs
            if "rank_001" in predicted_pdb.stem
        ]

    rmsds = []
    for predicted_pdb in predicted_pdbs:
        predicted = load_structure(predicted_pdb)
        predicted = predicted[np.isin(predicted.atom_name, "CA")]

        # iteratively superimpose each chain to get the best superimposition
        # TODO: will this work for a monomer?
        template, predicted = iteratively_superimpose_chains(template, predicted)
        rmsds.append(bts.rmsd(template, predicted))

    return float(np.mean(rmsds))


def compute_hydrophobicity(pdb_file_or_result_dir: Path) -> float:
    """
    Compute the fraction of the total SASA that is made up of hydrophobic atoms.
    """
    if pdb_file_or_result_dir.is_file():
        structure = load_structure(pdb_file_or_result_dir)
    else:
        top_pdb = [pdb for pdb in list(pdb_file_or_result_dir.glob("*.pdb")) if "rank_001" in pdb.stem][0]
        structure = load_structure(top_pdb)

    # consider only heavy atoms, otherwise we'd need to know what an H was bound to
    structure = structure[structure.element != "H"]

    sasa = bts.sasa(structure)

    elements = structure.element
    hydrophobic_sasa = [
        area if element in HYDROPHOBIC_ELEMENTS else 0
        for area, element in zip(sasa, elements)
    ]

    total_sasa = np.sum(sasa)
    total_hydrophobic_sasa = np.sum(hydrophobic_sasa)

    return total_hydrophobic_sasa / total_sasa


def make_chain_moved_pdb(
    original_pdb: Path, moved_pdb_directory: Path, distance: float = 100
) -> Path:
    """
    Take the original PDB, and move the first chain `distance` angstroms along the vector
    connecting the COM of whole PDB to the COM of the chain.
    """
    structure = load_structure(original_pdb)

    chains = bts.get_chains(structure)
    if len(chains) == 1:
        # don't try to move anything since this is just a monomer
        return original_pdb

    mobile_chain = chains[0]
    mobile_structure = structure[structure.chain_id == mobile_chain]
    fixed_structure = structure[structure.chain_id != mobile_chain]

    all_com = bts.centroid(structure)
    mobile_com = bts.centroid(mobile_structure)
    move_vector = mobile_com - all_com

    magnitude = np.linalg.norm(move_vector)
    if magnitude == 0:
        raise ValueError(
            f"{original_pdb} has a first chain whose COM overlaps with the full structure COM"
        )
    move_unit_vector = move_vector / magnitude

    moved_mobile = bts.translate(mobile_structure, move_unit_vector * distance)

    chain_moved_structure = fixed_structure + moved_mobile
    save_path = moved_pdb_directory / original_pdb.name
    save_structure(save_path, chain_moved_structure)

    return save_path
