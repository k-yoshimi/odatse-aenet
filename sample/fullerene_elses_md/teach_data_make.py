"""Extract trajectory data into XSF format for AENET training."""

import os

from ase.io import Trajectory


def extract_data(input_file, output_dir):
    """Convert an ASE trajectory to individual XSF files for AENET training.

    For each frame in the trajectory, an XSF file is written containing the
    total energy, cell vectors (for periodic systems), atomic positions, and
    forces. Periodic systems use ``CRYSTAL``/``PRIMVEC``/``PRIMCOORD`` format;
    non-periodic systems use ``ATOMS`` format.

    Parameters
    ----------
    input_file : str
        Path to the ASE trajectory file (``.traj``).
    output_dir : str
        Directory to write XSF files. Created if it does not exist.
        Files are named ``teach_data_1.xsf``, ``teach_data_2.xsf``, etc.
    """
    os.makedirs(output_dir, exist_ok=True)

    structures = Trajectory(input_file)

    for i, atoms in enumerate(structures):
        num_atoms = len(atoms)
        total_energy = atoms.get_total_energy()
        cell = atoms.get_cell()
        pbc = atoms.get_pbc()
        positions = atoms.get_positions()
        forces = atoms.get_forces()

        output_file = os.path.join(output_dir, f'teach_data_{i + 1}.xsf')

        with open(output_file, 'w') as f:
            if total_energy is not None:
                f.write(f"# total energy = {total_energy} eV\n\n")

            if all(pbc):
                f.write("CRYSTAL\n")
                f.write("PRIMVEC\n")
                for vec in cell:
                    f.write(f"  {vec[0]:.8f}  {vec[1]:.8f}  {vec[2]:.8f}\n")
                f.write("PRIMCOORD\n")
                f.write(f"{num_atoms} 1\n")
            elif not any(pbc):
                f.write("ATOMS\n")

            for k in range(num_atoms):
                symbol = atoms[k].symbol
                x, y, z = positions[k]
                fx, fy, fz = forces[k]
                f.write(f"{symbol} {x:.8f} {y:.8f} {z:.8f} {fx:.8f} {fy:.8f} {fz:.8f}\n")


if __name__ == "__main__":
    extract_data('benzene_optimization.traj', 'teach_data')
