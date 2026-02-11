"""Batch-run histogram plotting for all PAMC temperature step files."""

import argparse
import glob
import subprocess


def batch_plot(script, input_dir, title_option="tau", variable="x1"):
    """Run the histogram plotting script for each temperature step file.

    Parameters
    ----------
    script : str
        Path to ``plot_histogram.py``.
    input_dir : str
        Directory containing normalized weight summary files.
    title_option : str, optional
        Title label passed to the plotting script (default: ``'tau'``).
    variable : str, optional
        Variable to plot (default: ``'x1'``).
    """
    input_files = sorted(glob.glob(f"{input_dir}/*"))
    for i, input_file in enumerate(input_files):
        print(f"{i + 1}/{len(input_files)}: {input_file}")
        subprocess.run(["python3", script, input_file, str(i), title_option, variable])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch-run PAMC histogram plotting")
    parser.add_argument("script", help="Path to plot_histogram.py")
    parser.add_argument("input_dir", help="Directory with weight summary files")
    parser.add_argument("--title-option", default="tau", choices=["beta", "tau"],
                        help="Title label (default: tau)")
    parser.add_argument("--variable", default="x1",
                        choices=["x1", "x2", "x3", "both", "all"],
                        help="Variable(s) to plot (default: x1)")
    args = parser.parse_args()

    batch_plot(args.script, args.input_dir, args.title_option, args.variable)
