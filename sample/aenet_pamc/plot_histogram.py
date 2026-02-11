"""Plot 1D weighted histograms from PAMC output files."""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np


def plot_histograms(input_file, index, title_option, variables, bins=100, output_dir="1dhist"):
    """Plot weighted 1D histograms for selected variables.

    Parameters
    ----------
    input_file : str
        Path to the PAMC weight summary file.
    index : int
        Temperature step index, used for output filename.
    title_option : str
        Label for the title: ``'beta'`` or ``'tau'``.
    variables : list of str
        Variables to plot (e.g., ``['x1']``, ``['x1', 'x2', 'x3']``).
    bins : int, optional
        Number of histogram bins (default: 100).
    output_dir : str, optional
        Directory to save output images (default: ``'1dhist'``).
    """
    os.makedirs(output_dir, exist_ok=True)

    # Variable ranges: (min, max, tick_interval)
    var_ranges = {
        "x1": (0.2, 2.0, 0.1),
        "x2": (0.111281505036325, 0.888955610756928, 0.125),
        "x3": (0.110057883541601, 0.890182280028369, 0.125),
    }

    data = np.loadtxt(input_file, dtype="float", comments="#")
    weights = data[:, 4]
    weights = weights / np.sum(weights)

    # Read title value from the second line of the file
    with open(input_file, "r") as f:
        f.readline()
        second_line = f.readline()
    title_value = float(second_line.split()[3])

    fig, axes = plt.subplots(len(variables), 1, figsize=(8, len(variables) * 4))
    if len(variables) == 1:
        axes = [axes]

    for ax, var_name in zip(axes, variables):
        var_min, var_max, var_tick = var_ranges[var_name]
        var_index = list(var_ranges.keys()).index(var_name) + 5

        hist, bin_edges = np.histogram(
            data[:, var_index], bins=bins, range=(var_min, var_max), weights=weights
        )
        ax.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), edgecolor="black", align="edge")
        ax.set_xlabel(var_name)
        ax.set_ylabel("Weight")
        ax.set_title(f"1D Histogram of {var_name}_{title_option}_{title_value}")
        ax.set_xticks(np.arange(var_min, var_max + var_tick, var_tick))
        ax.grid()

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{index:03}_hist_{title_option}_{title_value}.png"))
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot 1D weighted histograms from PAMC output")
    parser.add_argument("input_file", help="Path to PAMC weight summary file")
    parser.add_argument("index", type=int, help="Temperature step index")
    parser.add_argument("title_option", choices=["beta", "tau"], help="Title label (beta or tau)")
    parser.add_argument(
        "variables",
        choices=["x1", "x2", "x3", "both", "all"],
        help="Variables to plot",
    )
    parser.add_argument("--bins", type=int, default=100, help="Number of histogram bins")
    parser.add_argument("--output-dir", default="1dhist", help="Output directory for plots")
    args = parser.parse_args()

    var_map = {
        "x1": ["x1"],
        "x2": ["x2"],
        "x3": ["x3"],
        "both": ["x1", "x2"],
        "all": ["x1", "x2", "x3"],
    }
    plot_histograms(
        args.input_file, args.index, args.title_option,
        var_map[args.variables], bins=args.bins, output_dir=args.output_dir,
    )
