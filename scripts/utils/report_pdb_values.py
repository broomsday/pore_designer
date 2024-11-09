"""
Load a PDB and report some design-related values.
"""

from pathlib import Path
from typing import Any

import typer
import pandas as pd
import numpy as np

import biotite.structure as bts

from designer.pdb import load_pdb, get_multimer_state, get_by_chain_sequence, compute_hydrophobicity, iteratively_superimpose_chains


def compute_rmsd(query: bts.AtomArray, template: bts.AtomArray) -> float:
    """
    Get the RMSD to the template PDB.
    """
    # load and reduce to Ca only
    template = template[np.isin(template.atom_name, "CA")]
    predicted = predicted[np.isin(predicted.atom_name, "CA")]

    template, predicted = iteratively_superimpose_chains(template, predicted)
    
    return bts.rmsd(template, predicted)


def get_plddt_from_structure(pdb: Path) -> float:
    """
    Compute the mean pLDDT from a structure where those values are encoded in the beta-field.
    """
    with open(pdb, mode="r", encoding="utf-8") as pdb_file:
        atom_lines = [line for line in pdb_file.readlines() if ("ATOM" in line) or ("HETATM" in line)]

    plddts = [float(line[61:67]) for line in atom_lines]

    return np.mean(plddts)


def analyze_pdb(pdb: Path, template_pdb: Path) -> dict[str, Any]:
    """
    Load a PDB file and report:

    1. the pLDDT
    2. the hydrophobicity
    3. the monomer length
    4. the number of chains
    5. (optional) RMSD to a template
    """
    query = load_pdb(pdb)

    plddt = get_plddt_from_structure(pdb)
    hydrophobicity = compute_hydrophobicity(pdb)
    sequence_length = len(get_by_chain_sequence(query)[0])
    multimer = get_multimer_state(query)

    if template_pdb is not None:
        template = load_pdb(template_pdb)
        rmsd = compute_rmsd(query, template)
    else:
        rmsd = None

    return {
        "pdb": pdb.stem,
        "plddt": f"{plddt:.2f}",
        "hydrophobicity": f"{hydrophobicity:.2f}",
        "monomer length": sequence_length,
        "oligomer state": multimer,
        "template_pdb": template_pdb,
        "rmsd": rmsd
    }


def main(
    query_pdb: Path = typer.Argument(..., help="PDB to analyze, or directory of PDBs."),
    template_pdb: Path = typer.Option(default=None, help="Template PDB to compute RMSD to if desired"),
):
    """
    For one or more PDBs report:

    1. the pLDDT
    2. the hydrophobicity
    3. the monomer length
    4. the number of chains
    5. (optional) RMSD to a template
    """
    if query_pdb.is_file():
        pdbs = [query_pdb]
    elif query_pdb.is_dir():
        pdbs = list(query_pdb.glob("*/*.pdb"))

    data_records = []
    for pdb in pdbs:
        data_records.append(analyze_pdb(pdb, template_pdb))

    data_df = pd.DataFrame().from_records(data_records)
    print(data_df)


if __name__ == "__main__":
    typer.run(main)

