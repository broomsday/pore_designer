"""
Determine pLDDT values across a number of oligomer states for a set of designs.

Additionally, assuming there are multiple seeds used, plot the pLDDT vs. oligomer state with uncertainty.
"""

from pathlib import Path

import typer
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from designer.file_utils import unzip_af2_prediction


def get_plddt(result: Path) -> float:
    """
    Pull out the mean pLDDT across the structure for one result.
    """
    result_df = pd.read_json(result)
    mean_plddt = np.mean(result_df["plddt"])

    return mean_plddt


def get_total_number_of_seeds(result_dir: Path) -> int:
    """
    Determine the number of unique seeds used.
    """
    seeds = set()
    for result in result_dir.glob("*scores*.json"):
        seed = result.stem.split("_")[-1]
        seeds.add(seed)

    return len(seeds)


def get_top_plddt_per_seed(result_dir: Path) -> list[float]:
    """
    Compute the mean pLDDT across the structure for the top prediction for each seed.
    """
    num_seeds = get_total_number_of_seeds(result_dir)

    seed_top_plddts = []
    for seed in range(num_seeds):
        seed_string = f"seed_{seed:03}"
        seed_plddts = []
        for result in result_dir.glob(f"*scores*{seed_string}.json"):
            seed_plddts.append(get_plddt(result))
        seed_top_plddts.append(max(seed_plddts))

    return seed_top_plddts


def plot_multiseed_plddts(data: pd.DataFrame, save: Path):
    """
    Plot each data-point as well as the mean.
    """
    plt.tight_layout()

    plt.scatter(data["oligomer"], data["plddt"], alpha=0.5)

    mean_values = data.groupby("oligomer")["plddt"].mean().reset_index()
    plt.scatter(mean_values["oligomer"], mean_values["plddt"], color="red")

    plt.title(data["name"].to_list()[0])
    plt.xlabel("Oligomer")
    plt.ylabel("pLDDT")
    plt.ylim(30, 100)

    plt.savefig(save)
    plt.close()


def main(
    input_dir: Path = typer.Argument(...),
    output_dir: Path = typer.Argument(...),
):
    """
    Perform AF2 noise analysis:

    1. Unzip all predictions in input_dir
    2. Pull out mean pLDDT from the top model for all predictions
    3. Plot mean pLDDT vs. oligomer state for all predictions with uncertainties

    Example:
    python scripts/analyze_af2_noise.py /home/broom/AlphaCarbon/code/pore_oligomer/data/design_r3/aron_shortlisted/final_aron_selections/af2_noise_test_predictions/ /home/broom/AlphaCarbon/code/pore_oligomer/data/design_r3/aron_shortlisted/final_aron_selections/af2_noise_test_analysis/
    """
    output_dir.mkdir(exist_ok=True)

    for archive in input_dir.glob("*.zip"):
        unzip_af2_prediction(archive, output_dir / archive.stem.replace(".result", ""), cleanup=False)

    data_dict = {"name": [], "oligomer": [], "plddt": []}
    for result_dir in output_dir.glob("*"):
        plddts = get_top_plddt_per_seed(result_dir)
        for plddt in plddts:
            data_dict["name"].append("_".join(result_dir.name.split("_")[:-1]))
            data_dict["oligomer"].append(int(result_dir.name.split("_")[-1]))
            data_dict["plddt"].append(plddt)
    data_df = pd.DataFrame().from_dict(data_dict)

    plot_dir = output_dir / "plots"
    plot_dir.mkdir(exist_ok=True)
    designs = {name for name in data_df.name}
    for design in designs:
        plot_multiseed_plddts(data_df[data_df.name == design], plot_dir / f"{design}_plddt.png")


if __name__ == "__main__":
    typer.run(main)
