"""Plot N-N distance vs total energy for QE reference and ANN predictions.

Compares three datasets:
1. ANN potential predictions (from ``predict.out``)
2. QE training data (from XSF files in relax directories)
3. QE reference data (from a separate XSF directory)
"""

import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np


def parse_xsf_files(directory):
    """Parse XSF files to extract N-N distances and total energies.

    Parameters
    ----------
    directory : str
        Directory containing ``.xsf`` files.

    Returns
    -------
    distances : list of float
        Inter-atomic distances in Angstrom.
    energies : list of float
        Total energies in eV.
    filenames : list of str
        Corresponding file paths.
    """
    distances = []
    energies = []
    filenames = []

    files = sorted(
        [f for f in os.listdir(directory) if f.endswith(".xsf")],
        key=lambda x: int(re.search(r"(\d+)", x).group(1)) if re.search(r"(\d+)", x) else 0,
    )

    for filename in files:
        filepath = os.path.join(directory, filename)
        with open(filepath, "r") as f:
            content = f.read()

        energy_match = re.search(r"total energy = ([\d.\-e]+) eV", content)
        if not energy_match:
            continue

        energy = float(energy_match.group(1))
        coordinates = re.findall(r"N\s+([\d.\-e]+)\s+([\d.\-e]+)\s+([\d.\-e]+)", content)

        if len(coordinates) >= 2:
            atom1 = np.array([float(x) for x in coordinates[0]])
            atom2 = np.array([float(x) for x in coordinates[1]])
            distance = np.linalg.norm(atom2 - atom1)

            distances.append(distance)
            energies.append(energy)
            filenames.append(filepath)

    return distances, energies, filenames


def parse_predict_out(filepath):
    """Parse AENET ``predict.out`` to extract distances and energies.

    Parameters
    ----------
    filepath : str
        Path to ``predict.out``.

    Returns
    -------
    distances : list of float
        Inter-atomic distances in Angstrom.
    energies : list of float
        Total energies in eV.
    """
    with open(filepath, "r") as f:
        content = f.read()

    blocks = content.split("+" * 70)
    distances = []
    energies = []

    for block in blocks:
        coordinates = re.findall(r"N\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", block)
        energy_match = re.search(r"Total energy\s+:\s+([-.\d]+) eV", block)

        if energy_match and len(coordinates) >= 2:
            atom1 = np.array([float(x) for x in coordinates[0]])
            atom2 = np.array([float(x) for x in coordinates[-1]])
            distance = np.linalg.norm(atom2 - atom1)

            distances.append(distance)
            energies.append(float(energy_match.group(1)))

    return distances, energies


def parse_test_set(train_out):
    """Extract test set indices from AENET ``train.out``.

    Parameters
    ----------
    train_out : str
        Path to ``train.out``.

    Returns
    -------
    indices : list of int
        Zero-based indices of test data points.
    """
    with open(train_out, "r", encoding="utf-8") as f:
        text = f.read()

    start = text.find(" Testing set : ")
    end = text.find("-" * 70, start)
    if start == -1 or end == -1:
        return []

    section = text[start:end]
    return [int(num) - 1 for num in re.findall(r"\d+", section)]


def plot_comparison(distances_ann, energies_ann,
                    distances_qe, energies_qe,
                    distances_test, energies_test,
                    output_file="distance_energy_plot.png",
                    xlim=(0.5, 2.0), ylim=None):
    """Plot distance-energy comparison between ANN and QE.

    Parameters
    ----------
    distances_ann : list of float
        ANN prediction distances.
    energies_ann : list of float
        ANN prediction energies.
    distances_qe : list of float
        QE training data distances.
    energies_qe : list of float
        QE training data energies.
    distances_test : list of float
        QE test data distances.
    energies_test : list of float
        QE test data energies.
    output_file : str, optional
        Output image path (default: ``'distance_energy_plot.png'``).
    xlim : tuple of float, optional
        X-axis limits (default: ``(0.5, 2.0)``).
    ylim : tuple of float, optional
        Y-axis limits (default: auto).
    """
    plt.figure(figsize=(10, 6))
    plt.scatter(distances_ann, energies_ann, color="b", label="ANN prediction", s=10)
    plt.scatter(distances_qe, energies_qe, color="r", label="QE training")
    if distances_test:
        plt.scatter(distances_test, energies_test, color="g", label="QE test", s=50)
    plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    plt.xlabel("Atomic Distance (Angstrom)")
    plt.ylabel("Total Energy (eV)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot N-N distance vs total energy (QE vs ANN)")
    parser.add_argument("--predict-out", default="predict.out",
                        help="Path to AENET predict.out")
    parser.add_argument("--train-out", default="../train/train.out",
                        help="Path to AENET train.out (for test set indices)")
    parser.add_argument("--qe-dir", default="../relax",
                        help="Base directory for QE training data")
    parser.add_argument("--n-structures", type=int, default=20,
                        help="Number of training structures")
    parser.add_argument("--output", default="distance_energy_plot.png",
                        help="Output image file")
    args = parser.parse_args()

    # ANN predictions
    distances_ann, energies_ann = parse_predict_out(args.predict_out)

    # QE training data from relax directories
    distances_qe = []
    energies_qe = []
    filenames_qe = []
    for i in range(args.n_structures):
        teach_dir = os.path.join(args.qe_dir, f"directory_{i}", "teach_data")
        if os.path.isdir(teach_dir):
            d, e, fn = parse_xsf_files(teach_dir)
            distances_qe.extend(d)
            energies_qe.extend(e)
            filenames_qe.extend(fn)

    # Test set
    distances_test = []
    energies_test = []
    if os.path.exists(args.train_out):
        test_indices = parse_test_set(args.train_out)
        for idx in test_indices:
            if idx < len(distances_qe):
                distances_test.append(distances_qe[idx])
                energies_test.append(energies_qe[idx])

    # Full range plot
    plot_comparison(
        distances_ann, energies_ann,
        distances_qe, energies_qe,
        distances_test, energies_test,
        output_file=args.output,
    )

    # Zoomed-in plot around the potential minimum
    base, ext = os.path.splitext(args.output)
    plot_comparison(
        distances_ann, energies_ann,
        distances_qe, energies_qe,
        distances_test, energies_test,
        output_file=f"{base}_zoom{ext}",
        ylim=(-770, -760),
    )
