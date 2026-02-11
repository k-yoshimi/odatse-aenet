"""
Unit tests for sample/fullerene_elses_md/teach_data_make.py

Tests cover:
1. extract_data function: XSF file generation from ASE Trajectory
2. Output format for periodic (CRYSTAL) and non-periodic (ATOMS) systems
3. Energy, positions, forces in output
"""

import os
import sys
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import Trajectory
from ase.calculators.singlepoint import SinglePointCalculator

sys.path.insert(0, str(Path(__file__).parent.parent / "sample" / "fullerene_elses_md"))

from teach_data_make import extract_data


def make_trajectory(tmpdir, atoms_list):
    """Create a trajectory file from list of (atoms, energy, forces) tuples."""
    traj_path = os.path.join(tmpdir, "test.traj")
    traj = Trajectory(traj_path, "w")
    for atoms, energy, forces in atoms_list:
        calc = SinglePointCalculator(atoms, energy=energy, forces=forces)
        atoms.calc = calc
        traj.write(atoms)
    traj.close()
    return traj_path


class TestExtractData(unittest.TestCase):
    """Test extract_data function."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_creates_output_directory(self):
        """Test that output directory is created."""
        atoms = Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.1]])
        forces = np.zeros((2, 3))
        traj_path = make_trajectory(self.tmpdir, [(atoms, -10.0, forces)])

        extract_data(traj_path, "output_xsf")
        self.assertTrue(os.path.isdir("output_xsf"))

    def test_creates_xsf_files(self):
        """Test that XSF files are created for each frame."""
        atoms = Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.1]])
        forces = np.zeros((2, 3))
        frames = [(atoms.copy(), -10.0, forces), (atoms.copy(), -11.0, forces)]
        traj_path = make_trajectory(self.tmpdir, frames)

        extract_data(traj_path, "output_xsf")
        self.assertTrue(os.path.exists("output_xsf/teach_data_1.xsf"))
        self.assertTrue(os.path.exists("output_xsf/teach_data_2.xsf"))

    def test_periodic_crystal_format(self):
        """Test output format for periodic systems (CRYSTAL + PRIMVEC)."""
        atoms = Atoms("N2",
                      positions=[[0, 0, 0], [0, 0, 1.1]],
                      cell=[10, 10, 10],
                      pbc=True)
        forces = np.array([[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]])
        traj_path = make_trajectory(self.tmpdir, [(atoms, -10.0, forces)])

        extract_data(traj_path, "output_xsf")

        with open("output_xsf/teach_data_1.xsf") as f:
            content = f.read()
        self.assertIn("CRYSTAL", content)
        self.assertIn("PRIMVEC", content)
        self.assertIn("PRIMCOORD", content)
        self.assertIn("2 1", content)  # natoms 1

    def test_non_periodic_atoms_format(self):
        """Test output format for non-periodic systems (ATOMS)."""
        atoms = Atoms("N2",
                      positions=[[0, 0, 0], [0, 0, 1.1]],
                      pbc=False)
        forces = np.zeros((2, 3))
        traj_path = make_trajectory(self.tmpdir, [(atoms, -10.0, forces)])

        extract_data(traj_path, "output_xsf")

        with open("output_xsf/teach_data_1.xsf") as f:
            content = f.read()
        self.assertIn("ATOMS", content)
        self.assertNotIn("CRYSTAL", content)

    def test_energy_in_header(self):
        """Test that total energy is written in the XSF header."""
        atoms = Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.1]], pbc=False)
        forces = np.zeros((2, 3))
        traj_path = make_trajectory(self.tmpdir, [(atoms, -12.345, forces)])

        extract_data(traj_path, "output_xsf")

        with open("output_xsf/teach_data_1.xsf") as f:
            content = f.read()
        self.assertIn("-12.345", content)
        self.assertIn("total energy", content)

    def test_positions_and_forces_written(self):
        """Test that positions and forces are correctly written."""
        pos = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        forces = [[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3]]
        atoms = Atoms("N2", positions=pos, pbc=False)
        traj_path = make_trajectory(self.tmpdir, [(atoms, -10.0, np.array(forces))])

        extract_data(traj_path, "output_xsf")

        with open("output_xsf/teach_data_1.xsf") as f:
            lines = f.readlines()
        # Find atom lines (after ATOMS line)
        atom_lines = [l for l in lines if l.startswith("N ")]
        self.assertEqual(len(atom_lines), 2)

        parts = atom_lines[0].split()
        self.assertAlmostEqual(float(parts[1]), 1.0, places=5)
        self.assertAlmostEqual(float(parts[4]), 0.1, places=5)

    def test_existing_output_dir(self):
        """Test that function works when output directory already exists."""
        os.makedirs("output_xsf")
        atoms = Atoms("N2", positions=[[0, 0, 0], [0, 0, 1.1]], pbc=False)
        forces = np.zeros((2, 3))
        traj_path = make_trajectory(self.tmpdir, [(atoms, -10.0, forces)])

        # Should not raise
        extract_data(traj_path, "output_xsf")
        self.assertTrue(os.path.exists("output_xsf/teach_data_1.xsf"))


if __name__ == "__main__":
    unittest.main()
