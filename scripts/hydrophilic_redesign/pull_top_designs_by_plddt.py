"""
Pull out the top N designs by pLDDT.
"""

from pathlib import Path
import shutil

import typer
import numpy as np
import pandas as pd
from tqdm import tqdm


MIN_PLDDT = 85


def get_plddt(result: Path) -> float:
    """
    Pull out the mean pLDDT across the structure for one result.
    """
    result_df = pd.read_json(result)
    mean_plddt = np.mean(result_df["plddt"])

    return mean_plddt


def main(
    af2_result_dir: Path = typer.Argument(..., help="The directory containing sub-dirs of results."),
    out_dir: Path = typer.Argument(..., help="Where to save the winners"),
    min_plddt: float = typer.Option(default=MIN_PLDDT),
    dry_run: bool = typer.Option(default=True),
):
    """
    Look through the AF design results in `design_dir` and pull out the top `num` designs by pLDDT.
     
    Dump the structures and a .csv of structure,plddt in the `out_dir`
    """
    design_dict = {"id": [], "plddt": [], "pdb": []}
    for prediction_dir in tqdm(list(af2_result_dir.glob("*")), desc="Pulling data"):
        if not prediction_dir.is_dir():
            continue

        top_json = list(prediction_dir.glob("*rank_001*.json"))[0]
        plddt = get_plddt(top_json)

        if plddt > min_plddt:
            top_pdb = list(prediction_dir.glob("*rank_001*.pdb"))[0]

            design_dict["id"].append(prediction_dir.name)
            design_dict["plddt"].append(plddt)
            design_dict["pdb"].append(top_pdb.name)

    design_df = pd.DataFrame().from_dict(design_dict)
    design_df = design_df.sort_values(by="plddt", ascending=False)
    design_df = design_df[design_df["plddt"] >= min_plddt]

    out_dir.mkdir(exist_ok=True, parents=True)
    for _, design_row in design_df.iterrows():
        source_pdb = af2_result_dir / design_row.id / design_row.pdb
        destintation_pdb = out_dir / f"{design_row.id}.pdb"
        if not dry_run:
            shutil.copy(source_pdb, destintation_pdb)

    design_df.to_csv(out_dir / "summary.csv")
    print(design_df)
    print(f"Found {len(design_df.index)} high pLDDT designs")
    if dry_run:
        print(f"Nothing copied as this is a dry-run")


if __name__ == "__main__":
    typer.run(main)