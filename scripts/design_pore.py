"""
Read a 'pore designer' config file to launch of continue a pore design job.
"""

from pathlib import Path
import json

import typer
import pandas as pd
from biotite.structure.io import save_structure

from designer import design, pdb


app = typer.Typer()


@app.command()
def make_config_positive(
    output_dir: Path = typer.Argument(
        ..., help="Directory to save config and job output"
    ),
    input_pdb: Path = typer.Argument(..., help="PDB to design into a pore."),
    num_mpnn: int = typer.Option(
        default=10000, help="Number of ProteinMPNN designs to make"
    ),
    temperature_mpnn: float = typer.Option(
        default=0.1, help="MPNN sampling temperature."
    ),
    num_af2: int = typer.Option(
        default=100, help="Number of design sequences to test by AlphaFold2"
    ),
    recycle_af2: int = typer.Option(default=3, help="Number of recycles to use in AF2"),
    top_plddt: float = typer.Option(
        default=90, help="AF2 pLDDT cutoff to select final sequences"
    ),
    mean_plddt: float = typer.Option(
        default=80, help="AF2 pLDDT cutoff to select final sequences"
    ),
    top_rmsd: float = typer.Option(default=1.0, help="Max RMSD for selection"),
    mean_rmsd: float = typer.Option(default=2.0, help="Max RMSD for selection"),
    select_identity: float = typer.Option(
        default=0.9, help="Maximum identity between any two selected sequences"
    ),
    select_oligomer_rank: int = typer.Option(
        default=1,
        help="Minimum oligomer check rank to be selected for oligomer designs",
    ),
    oligomer_lower_offset: int = typer.Option(
        default=3,
        help="How many fewer oligomers to test when performing the oligomer check",
    ),
    oligomer_higher_offset: int = typer.Option(
        default=4,
        help="How many more oligomers to test when performing the oligomer check",
    ),
    symmetry_dict: Path = typer.Option(
        default=None,
        help="If using a monomer from 'monomerizer' provide the symmetry dict output",
    ),
    fixed_dict: Path = typer.Option(
        default=None,
        help="If using a monomer from 'monomerizer' optionally provide a fixed dict output",
    ),
    overwrite: bool = typer.Option(default=False, help="Overwrite existing config"),
):
    """
    Based on CLI input, generate a config file for a given positive design project.
    """
    if output_dir.is_dir():
        if not overwrite:
            raise FileExistsError("Config already  exists, use --overwrite")
        else:
            print("Directory already exists, clean manually if you want a fresh run")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "input").mkdir(exist_ok=True)

    # clean the input structure
    clean_input_structure = pdb.clean_pdb(input_pdb)
    clean_input_pdb = output_dir / "input" / "input.pdb"
    save_structure(clean_input_pdb, clean_input_structure)

    values = {
        "job_type": "positive",
        "directory": str(output_dir.absolute()),
        "input_pdb": str(clean_input_pdb.absolute()),
        "num_mpnn": num_mpnn,
        "temperature_mpnn": temperature_mpnn,
        "num_af2": num_af2,
        "recycle_af2": recycle_af2,
        "top_plddt": top_plddt,
        "mean_plddt": mean_plddt,
        "top_rmsd": top_rmsd,
        "mean_rmsd": mean_rmsd,
        "select_identity": select_identity,
        "select_oligomer_rank": select_oligomer_rank,
        "oligomer_lower_offset": oligomer_lower_offset,
        "oligomer_higher_offset": oligomer_higher_offset,
        "multimer": pdb.get_multimer_state(clean_input_structure),
    }

    if symmetry_dict is None:
        values["symmetry_dict"] = None
    else:
        values["symmetry_dict"] = str(symmetry_dict.absolute())

    if fixed_dict is None:
        values["fixed_dict"] = None
    else:
        values["fixed_dict"] = str(fixed_dict.absolute())

    with open(output_dir / "config.json", mode="w", encoding="utf-8") as config_file:
        json.dump(values, config_file)


@app.command()
def make_config_negative(
    output_dir: Path = typer.Argument(
        ..., help="Directory to save config and job output"
    ),
    positive_pdb: Path = typer.Argument(..., help="Positive design PDB."),
    negative_pdb_dir: Path = typer.Argument(
        ..., help="Directory of all negative PDBs."
    ),
    num_mpnn: int = typer.Option(
        default=10000, help="Number of ProteinMPNN designs to make"
    ),
    temperature_mpnn: float = typer.Option(
        default=0.1, help="MPNN sampling temperature."
    ),
    num_af2: int = typer.Option(
        default=100, help="Number of design sequences to test by AlphaFold2"
    ),
    recycle_af2: int = typer.Option(default=3, help="Number of recycles to use in AF2"),
    top_plddt: float = typer.Option(
        default=90, help="AF2 pLDDT cutoff to select final sequences"
    ),
    mean_plddt: float = typer.Option(
        default=80, help="AF2 pLDDT cutoff to select final sequences"
    ),
    top_rmsd: float = typer.Option(default=1.0, help="Max RMSD for selection"),
    mean_rmsd: float = typer.Option(default=2.0, help="Max RMSD for selection"),
    select_identity: float = typer.Option(
        default=0.9, help="Maximum identity between any two selected sequences"
    ),
    select_oligomer_rank: int = typer.Option(
        default=1,
        help="Minimum oligomer check rank to be selected for oligomer designs",
    ),
    oligomer_lower_offset: int = typer.Option(
        default=3,
        help="How many fewer oligomers to test when performing the oligomer check",
    ),
    oligomer_higher_offset: int = typer.Option(
        default=4,
        help="How many more oligomers to test when performing the oligomer check",
    ),
    overwrite: bool = typer.Option(default=False, help="Overwrite existing config"),
):
    """
    Based on CLI input, generate a config file for a given negative design project.
    """
    if output_dir.is_dir():
        if not overwrite:
            raise FileExistsError("Config already  exists, use --overwrite")
        else:
            print("Directory already exists, clean manually if you want a fresh run")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "input").mkdir(exist_ok=True)
    (output_dir / "input" / "positive").mkdir(exist_ok=True)
    (output_dir / "input" / "negative").mkdir(exist_ok=True)

    # clean the positive input structure
    clean_positive_structure = pdb.clean_pdb(positive_pdb)
    clean_positive_pdb = output_dir / "input" / "positive" / positive_pdb.name
    save_structure(clean_positive_pdb, clean_positive_structure)

    # clean the negative input structures
    clean_negative_pdbs = output_dir / "input" / "negative"
    for negative_pdb in negative_pdb_dir.glob("*.pdb"):
        clean_negative_structure = pdb.clean_pdb(negative_pdb)
        clean_negative_pdb = output_dir / "input" / "negative" / negative_pdb.name
        save_structure(clean_negative_pdb, clean_negative_structure)

    values = {
        "job_type": "negative",
        "directory": str(output_dir.absolute()),
        "positive_pdb": str(clean_positive_pdb.absolute()),
        "negative_pdbs": str(clean_negative_pdbs.absolute()),
        "num_mpnn": num_mpnn,
        "temperature_mpnn": temperature_mpnn,
        "num_af2": num_af2,
        "recycle_af2": recycle_af2,
        "top_plddt": top_plddt,
        "mean_plddt": mean_plddt,
        "top_rmsd": top_rmsd,
        "mean_rmsd": mean_rmsd,
        "select_identity": select_identity,
        "select_oligomer_rank": select_oligomer_rank,
        "oligomer_lower_offset": oligomer_lower_offset,
        "oligomer_higher_offset": oligomer_higher_offset,
        "multimer": pdb.get_multimer_state(clean_positive_structure),
    }

    with open(output_dir / "config.json", mode="w", encoding="utf-8") as config_file:
        json.dump(values, config_file)


@app.command()
def make_config_mutation(
    output_dir: Path = typer.Argument(
        ..., help="Directory to save config and job output"
    ),
    input_pdb: Path = typer.Argument(..., help="PDB to design into a pore."),
    mutations_per_seq: int = typer.Option(
        default=1, help="How many mutations to make to each sequence"
    ),
    mutation_distribution: Path = typer.Option(
        default=None, help="Path to a .json distribution to pull mutations from"
    ),
    num_af2: int = typer.Option(
        default=100,
        help="Number of design sequences to test by AlphaFold2, same as the number of mutated sequences",
    ),
    recycle_af2: int = typer.Option(default=3, help="Number of recycles to use in AF2"),
    top_plddt: float = typer.Option(
        default=90, help="AF2 pLDDT cutoff to select final sequences"
    ),
    mean_plddt: float = typer.Option(
        default=80, help="AF2 pLDDT cutoff to select final sequences"
    ),
    top_rmsd: float = typer.Option(default=1.0, help="Max RMSD for selection"),
    mean_rmsd: float = typer.Option(default=2.0, help="Max RMSD for selection"),
    select_identity: float = typer.Option(
        default=0.9, help="Maximum identity between any two selected sequences"
    ),
    select_oligomer_rank: int = typer.Option(
        default=1,
        help="Minimum oligomer check rank to be selected for oligomer designs",
    ),
    oligomer_lower_offset: int = typer.Option(
        default=3,
        help="How many fewer oligomers to test when performing the oligomer check",
    ),
    oligomer_higher_offset: int = typer.Option(
        default=4,
        help="How many more oligomers to test when performing the oligomer check",
    ),
    overwrite: bool = typer.Option(default=False, help="Overwrite existing config"),
):
    """
    Based on CLI input, generate a config file for a given positive design project.
    """
    if output_dir.is_dir():
        if not overwrite:
            raise FileExistsError("Config already  exists, use --overwrite")
        else:
            print("Directory already exists, clean manually if you want a fresh run")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "input").mkdir(exist_ok=True)

    # clean the input structure
    clean_input_structure = pdb.clean_pdb(input_pdb)
    clean_input_pdb = output_dir / "input" / "input.pdb"
    save_structure(clean_input_pdb, clean_input_structure)

    values = {
        "job_type": "mutation",
        "directory": str(output_dir.absolute()),
        "input_pdb": str(clean_input_pdb.absolute()),
        "mutations_per_seq": mutations_per_seq,
        "mutation_distribution": str(mutation_distribution.absolute()),
        "num_af2": num_af2,
        "recycle_af2": recycle_af2,
        "top_plddt": top_plddt,
        "mean_plddt": mean_plddt,
        "top_rmsd": top_rmsd,
        "mean_rmsd": mean_rmsd,
        "select_identity": select_identity,
        "select_oligomer_rank": select_oligomer_rank,
        "oligomer_lower_offset": oligomer_lower_offset,
        "oligomer_higher_offset": oligomer_higher_offset,
        "multimer": pdb.get_multimer_state(clean_input_structure),
    }

    with open(output_dir / "config.json", mode="w", encoding="utf-8") as config_file:
        json.dump(values, config_file)


@app.command()
def make_config_metric(
    output_dir: Path = typer.Argument(
        ..., help="Directory to save config and job output"
    ),
    input_csv: Path = typer.Argument(..., help="Positive design PDB."),
    recycle_af2: int = typer.Option(default=3, help="Number of recycles to use in AF2"),
    oligomer_lower_offset: int = typer.Option(
        default=3,
        help="How many fewer oligomers to test when performing the oligomer check",
    ),
    oligomer_higher_offset: int = typer.Option(
        default=4,
        help="How many more oligomers to test when performing the oligomer check",
    ),
    overwrite: bool = typer.Option(default=False, help="Overwrite existing config"),
):
    """
    Based on CLI input, generate a config file for a given metric project.
    """
    if output_dir.is_dir():
        if not overwrite:
            raise FileExistsError("Config already  exists, use --overwrite")
        else:
            print("Directory already exists, clean manually if you want a fresh run")

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "input").mkdir(exist_ok=True)

    # copy over the input .csv
    input_data = pd.read_csv(input_csv, index_col="pdb")
    input_data.to_csv(output_dir / "input" / "input.csv")

    values = {
        "job_type": "metric",
        "directory": str(output_dir.absolute()),
        "input_csv": str((output_dir / "input" / "input.csv").absolute()),
        "recycle_af2": recycle_af2,
        "oligomer_lower_offset": oligomer_lower_offset,
        "oligomer_higher_offset": oligomer_higher_offset,
        "multimer": max(input_data["oligomer"]),
    }

    with open(output_dir / "config.json", mode="w", encoding="utf-8") as config_file:
        json.dump(values, config_file)


@app.command()
def design_pore(
    config_path: Path = typer.Argument(..., help="Pore Design config file")
):
    """
    Parse a `pore designer` config file and launch or continue the pore design job.
    """
    with open(config_path, mode="r", encoding="utf-8") as config_file:
        config = json.load(config_file)

    if config["job_type"] == "positive":
        design.design_pore_positive(config_path)
    elif config["job_type"] == "negative":
        design.design_pore_negative(config_path)
    elif config["job_type"] == "mutation":
        design.design_pore_mutations(config_path)
    elif config["job_type"] == "metric":
        design.produce_multimer_metrics(config_path)


if __name__ == "__main__":
    app()
