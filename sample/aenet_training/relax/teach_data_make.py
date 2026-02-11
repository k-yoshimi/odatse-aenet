"""Convert Quantum ESPRESSO output to XSF files for AENET training.

Parses a QE ``.pwo`` output file and extracts atomic positions, forces,
and total energies at each ionic step, writing them as individual XSF files.
"""

import argparse
import os
import re

from ase.units import Ry


def extract_data(input_file, output_dir):
    """Extract training data from a Quantum ESPRESSO output file.

    Parameters
    ----------
    input_file : str
        Path to the QE ``.pwo`` output file.
    output_dir : str
        Directory to write XSF files. Created if it does not exist.
        Files are named ``teach_data_1.xsf``, ``teach_data_2.xsf``, etc.

    Notes
    -----
    Energy is converted from Ry to eV using ``ase.units.Ry``.
    Forces are converted from Ry/Bohr to eV/Angstrom using the factor
    ``Ry / Bohr = 25.7110`` eV/Angstrom.
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(input_file, "r") as f:
        lines = f.readlines()

    total_energy_pattern = re.compile(r"!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry")
    atomic_positions_pattern = re.compile(r"ATOMIC_POSITIONS \(angstrom\)")
    force_pattern = re.compile(r"force\s+=\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)")

    # Ry/Bohr -> eV/Angstrom
    force_conversion = Ry / 0.529177  # ≈ 25.7110

    total_energy_list = []
    atomic_positions = []
    forces = []
    current_positions = []

    for line in lines:
        energy_match = total_energy_pattern.search(line)
        force_match = force_pattern.search(line)

        if energy_match:
            total_energy_list.append(Ry * float(energy_match.group(1)))

        elif atomic_positions_pattern.search(line):
            if current_positions:
                atomic_positions.append(current_positions)
                current_positions = []
            current_positions.append("ATOMS")

        elif current_positions:
            if line.strip() == "End final coordinates":
                atomic_positions.append(current_positions)
                current_positions = []
            elif "0   0   0" in line.strip():
                position = line.strip().replace("0   0   0", "").strip()
                current_positions.append(position)
            else:
                stripped = line.strip()
                if stripped:
                    current_positions.append(stripped)

        elif force_match:
            fx = force_conversion * float(force_match.group(1))
            fy = force_conversion * float(force_match.group(2))
            fz = force_conversion * float(force_match.group(3))
            forces.append(f"{fx:.8f} {fy:.8f} {fz:.8f}")

    if current_positions:
        atomic_positions.append(current_positions)

    # Number of atoms per frame (from first frame, excluding "ATOMS" header)
    natoms = len(atomic_positions[0]) - 1 if atomic_positions else 0

    for idx, positions in enumerate(atomic_positions):
        output_file = os.path.join(output_dir, f"teach_data_{idx + 1}.xsf")
        with open(output_file, "w") as f:
            if idx < len(total_energy_list):
                f.write(f"# total energy = {total_energy_list[idx]} eV\n\n")
            for i, position in enumerate(positions):
                if i == 0:
                    f.write(f"{position}\n")
                else:
                    force_idx = idx * natoms + (i - 1)
                    if force_idx < len(forces):
                        f.write(f"{position}     {forces[force_idx]}\n")
                    else:
                        f.write(f"{position}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert QE output to XSF files for AENET training"
    )
    parser.add_argument("--input", default="n2_dimer.pwo", help="QE output file (.pwo)")
    parser.add_argument("--output-dir", default="teach_data", help="Output directory for XSF files")
    args = parser.parse_args()

    extract_data(args.input, args.output_dir)
