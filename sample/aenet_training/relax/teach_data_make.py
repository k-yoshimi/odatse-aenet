"""Convert Quantum ESPRESSO output to XSF files for AENET training.

Reads a QE ``.pwo`` output file using ASE and extracts atomic positions,
forces, and total energies, writing them as individual XSF files.
"""

import argparse
import os

from ase.io import read


def extract_data(input_file, output_dir):
    """Extract training data from a Quantum ESPRESSO output file.

    Parameters
    ----------
    input_file : str
        Path to the QE ``.pwo`` output file.
    output_dir : str
        Directory to write XSF files. Created if it does not exist.
        Files are named ``teach_data_1.xsf``, ``teach_data_2.xsf``, etc.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Read all frames from QE output (SCF: 1 frame, relax: multiple frames)
    frames = read(input_file, index=":", format="espresso-out")
    if not isinstance(frames, list):
        frames = [frames]

    for idx, atoms in enumerate(frames):
        energy = atoms.get_potential_energy()  # eV
        forces = atoms.get_forces()  # eV/Angstrom
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()

        output_file = os.path.join(output_dir, f"teach_data_{idx + 1}.xsf")
        with open(output_file, "w") as f:
            f.write(f"# total energy = {energy} eV\n\n")
            f.write("ATOMS\n")
            for symbol, pos, frc in zip(symbols, positions, forces):
                f.write(
                    f"{symbol}  {pos[0]:.10f}  {pos[1]:.10f}  {pos[2]:.10f}"
                    f"     {frc[0]:.8f}  {frc[1]:.8f}  {frc[2]:.8f}\n"
                )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert QE output to XSF files for AENET training"
    )
    parser.add_argument("--input", default="n2_dimer.pwo", help="QE output file (.pwo)")
    parser.add_argument("--output-dir", default="teach_data", help="Output directory for XSF files")
    args = parser.parse_args()

    extract_data(args.input, args.output_dir)
