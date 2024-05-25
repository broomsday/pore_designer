from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from designer import alphafold


METRICS = {
    "top_plddt": "higher",
    "mean_plddt": "higher",
    "top_pae": "lower",
    "mean_pae": "lower",
    "top_max_pae": "lower",
    "mean_max_pae": "lower",
    "top_ptm": "higher",
    "mean_ptm": "higher",
    "top_iptm": "higher",
    "mean_iptm": "higher",
    "top_quad": "higher",
    "mean_quad": "higher",
    "top_mpnn": "lower",
}

FONTSIZE = 12
HIGHLIGHT_COLOR = "red"


# TODO: rename this if going beyond just pLDDT
def make_oligomer_v_plddt_plot(
    oligomer_values: pd.DataFrame,
    design_id: str,
    save_dir: Path,
    intended_oligomer: int,
    min_plddt: float = 30,
) -> None:
    """
    Produce a split plot of oligomer vs. plddt, one for top model only and one for all models.
    """
    # TODO: do 5 plots: top_plddt, top_pae, top_ptm, top_iptm, top_mpnn
    figure, axes = plt.subplots(2, sharex=True, sharey=True)
    axes[0].set_ylim(min_plddt, 100)
    figure.set_figwidth(5)
    figure.set_figheight(8)

    # plot top pLDDT and best oligomer
    axes[0].scatter(
        oligomer_values.oligomer.to_list(), oligomer_values.top_plddt.to_list()
    )
    best_oligomer, best_plddt = alphafold.get_best_oligomer(
        oligomer_values, "top_plddt"
    )
    axes[0].scatter(best_oligomer, best_plddt, c="red")

    # plot mean pLDDT and best oligomer
    axes[1].scatter(
        oligomer_values.oligomer.to_list(), oligomer_values.mean_plddt.to_list()
    )
    best_oligomer, best_plddt = alphafold.get_best_oligomer(
        oligomer_values, "mean_plddt"
    )
    axes[1].scatter(best_oligomer, best_plddt, c="red")

    if intended_oligomer is not None:
        axes[0].axvline(x=intended_oligomer, color="red")
        axes[1].axvline(x=intended_oligomer, color="red")

    axes[0].set_title(design_id)
    axes[0].set(ylabel="top model pLDDT")
    axes[1].set(xlabel="oligomers", ylabel="all models pLDDT")
    plt.tight_layout()
    plt.savefig(save_dir / f"{design_id}.png")
    plt.close()


def make_metric_correlation_plot(
    experimental_and_predicted: pd.DataFrame,
    metric: str,
    save_dir: Path,
) -> None:
    """
    Make a scatter plot with line of best fit and equation.
    """
    x = experimental_and_predicted["experimental"]
    y = experimental_and_predicted["predicted"]

    plt.scatter(x, y)

    # compute the line of best fit and plot
    slope, intercept = np.polyfit(x, y, 1)
    correlation = np.corrcoef(x, y)[0, 1]
    plt.plot(
        x,
        slope * x + intercept,
        color="red",
    )
    equation = f"y = {slope:.2f}x + {intercept:.2f}, r = {correlation:.2f}"
    plt.text(np.min(x), np.max(y), equation, fontsize=FONTSIZE, color=HIGHLIGHT_COLOR)

    plt.title(f"{metric} correlation")
    plt.xlabel("Experimental")
    plt.ylabel("Predicted")
    plt.tight_layout()
    plt.savefig(save_dir / f"{metric}.png")
    plt.close()


def plot_oligomer_check(config: dict) -> None:
    """
    Load the plddt values for each oligomer and plot.
    """
    # load the full and selected oligomer results
    oligomer_values = pd.read_csv(
        Path(config["directory"]) / "oligomer_values.csv", index_col=0
    )
    selected_seqs = alphafold.load_alphafold(config, "oligomer", "selected")

    # for each selected design, make the plots
    final_directory = Path(config["directory"]) / "final_selected"
    for selected in selected_seqs:
        seq_directory = final_directory / str(selected.id)
        selected_oligomer_values = oligomer_values[
            oligomer_values["design_id"] == selected.id
        ]

        make_oligomer_v_plddt_plot(
            selected_oligomer_values, selected.id, seq_directory, config["multimer"]
        )


def plot_all_oligomer_checks(config: dict) -> None:
    """
    Make an oligomer plot for all designs/PDBs.
    """
    # load the full and selected oligomer results
    oligomer_values = pd.read_csv(
        Path(config["directory"]) / "oligomer_values.csv", index_col=0
    )

    # set a column to hold the base design/PDB value without the oligomer notation
    oligomer_values["base_pdb"] = oligomer_values["design_id"].apply(
        lambda id: "_".join(id.split("_")[:-1])
    )
    base_pdbs = set(oligomer_values["base_pdb"].to_list())

    # iterate over all results for each base pdb
    final_directory = Path(config["directory"]) / "final_selected"
    final_directory.mkdir(exist_ok=True)
    for base_pdb in base_pdbs:
        #  pull out values for just this PDB
        selected_oligomer_values = oligomer_values[
            oligomer_values["base_pdb"] == base_pdb
        ]

        # use the WT column to get the WT oligomer value
        wt_oligomer = int(
            selected_oligomer_values[selected_oligomer_values["wt"]]["oligomer"]
        )

        # generate the plot
        seq_directory = final_directory / str(base_pdb)
        seq_directory.mkdir(exist_ok=True)
        make_oligomer_v_plddt_plot(
            selected_oligomer_values, base_pdb, seq_directory, wt_oligomer
        )


def get_top_predicted_by_metric(
    pdb: str, metric: str, oligomer_values: pd.DataFrame
) -> int:
    """
    For a given `pdb` and `metric`, find the oligomer value with the best metric.
    """
    pdb_values = oligomer_values[oligomer_values["pdb"] == pdb][
        ["oligomer", metric]
    ].sort_values(by=metric)

    if METRICS[metric] == "higher":
        return int(pdb_values.iloc[-1]["oligomer"])
    else:
        return int(pdb_values.iloc[0]["oligomer"])


def compute_metric_correlations(oligomer_values: pd.DataFrame) -> dict[str, float]:
    """
    Compute the correlation coefficient for experimental oligomer vs. predicted,
    based on each of the top metrics.
    """
    # add the PDB as a non-unique element for downstream grouping
    oligomer_values["pdb"] = oligomer_values["design_id"].apply(
        lambda id: id.split("_")[0]
    )
    # get the experimental oligomer values
    experimental_and_predicted = oligomer_values[oligomer_values["wt"]][
        ["pdb", "oligomer"]
    ]
    experimental_and_predicted = experimental_and_predicted.rename(
        columns={"oligomer": "experimental"}
    )

    metric_correlations = dict()
    for metric in METRICS:
        experimental_and_predicted["predicted"] = experimental_and_predicted[
            "pdb"
        ].apply(get_top_predicted_by_metric, args=([metric, oligomer_values]))
        metric_correlations[metric] = np.corrcoef(
            experimental_and_predicted["experimental"],
            experimental_and_predicted["predicted"],
        )[0, 1]

    return metric_correlations


def plot_metric_correlations(oligomer_values: pd.DataFrame, save_dir: Path) -> None:
    """
    Plot the true WT oligomer vs. the predicted based on each of the top metrics.
    """
    save_dir.mkdir(exist_ok=True, parents=True)
    # add the PDB as a non-unique element for downstream grouping
    oligomer_values["pdb"] = oligomer_values["design_id"].apply(
        lambda id: id.split("_")[0]
    )
    # get the experimental oligomer values
    experimental_and_predicted = oligomer_values[oligomer_values["wt"]][
        ["pdb", "oligomer"]
    ]
    experimental_and_predicted = experimental_and_predicted.rename(
        columns={"oligomer": "experimental"}
    )

    for metric in METRICS:
        experimental_and_predicted["predicted"] = experimental_and_predicted[
            "pdb"
        ].apply(get_top_predicted_by_metric, args=([metric, oligomer_values]))
        make_metric_correlation_plot(experimental_and_predicted, metric, save_dir)
