"""
Unit tests for ElsesCalculator.

Tests cover:
1. Module import and class hierarchy
2. Unit conversion constants
3. write_input: XYZ file generation from ASE Atoms
4. read_results: energy/force extraction from restart.xml and Output.txt
"""

import os
import sys
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.calculator import PropertyNotPresent

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ElsesCalculator.elses import elses, ev, angst, eV_angst

# Dummy command to satisfy ASE FileIOCalculator requirement
DUMMY_COMMAND = "echo dummy"


class TestImport(unittest.TestCase):
    """Test module import and class structure."""

    def test_import_elses(self):
        from ElsesCalculator import elses as cls
        self.assertIsNotNone(cls)

    def test_inherits_fileio_calculator(self):
        from ase.calculators.calculator import FileIOCalculator
        self.assertTrue(issubclass(elses, FileIOCalculator))

    def test_implemented_properties(self):
        self.assertIn("energy", elses.implemented_properties)
        self.assertIn("forces", elses.implemented_properties)


class TestUnitConstants(unittest.TestCase):
    """Test unit conversion constants."""

    def test_ev_hartree_to_ev(self):
        # 1 Hartree ≈ 27.2114 eV (from ase.units)
        from ase.units import Hartree
        self.assertAlmostEqual(ev, Hartree, places=8)

    def test_angst_bohr_to_angstrom(self):
        # 1 Bohr ≈ 0.52918 Angstrom (from ase.units)
        from ase.units import Bohr
        self.assertAlmostEqual(angst, Bohr, places=8)

    def test_eV_angst(self):
        self.assertAlmostEqual(eV_angst, ev / angst, places=5)


class TestWriteInput(unittest.TestCase):
    """Test write_input XYZ file generation."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_write_input_n2(self):
        """Test XYZ output for N2 molecule."""
        atoms = Atoms("N2",
                      positions=[[0, 0, 0], [0, 0, 1.1]],
                      cell=[10, 10, 10],
                      pbc=True)
        calc = elses(command=DUMMY_COMMAND)
        calc.write_input(atoms)

        self.assertTrue(os.path.exists("N2.xyz"))
        with open("N2.xyz") as f:
            lines = f.readlines()
        # First line: num_atoms + cell params
        self.assertTrue(lines[0].startswith("2"))
        self.assertIn("10.0", lines[0])
        # Second line: formula
        self.assertEqual(lines[1].strip(), "N2")
        # Atom lines
        self.assertEqual(len(lines), 4)  # header + formula + 2 atoms
        self.assertTrue(lines[2].strip().startswith("N"))
        self.assertTrue(lines[3].strip().startswith("N"))

    def test_write_input_single_element(self):
        """Test XYZ output for single-element case (C atom)."""
        atoms = Atoms("C",
                      positions=[[0, 0, 0]],
                      cell=[5, 5, 5],
                      pbc=True)
        calc = elses(command=DUMMY_COMMAND)
        calc.write_input(atoms)

        self.assertTrue(os.path.exists("C.xyz"))
        with open("C.xyz") as f:
            lines = f.readlines()
        self.assertTrue(lines[0].startswith("1"))
        self.assertEqual(lines[1].strip(), "C")

    def test_write_input_multi_species(self):
        """Test XYZ output for multi-species system (NaCl)."""
        atoms = Atoms("NaCl",
                      positions=[[0, 0, 0], [2.82, 0, 0]],
                      cell=[5.64, 5.64, 5.64],
                      pbc=True)
        calc = elses(command=DUMMY_COMMAND)
        calc.write_input(atoms)

        # Filename should contain both elements
        found = False
        for name in ["NaCl.xyz", "ClNa.xyz"]:
            if os.path.exists(name):
                found = True
                break
        self.assertTrue(found)

    def test_write_input_positions(self):
        """Test that positions are written correctly."""
        pos = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        atoms = Atoms("N2",
                      positions=pos,
                      cell=[10, 10, 10],
                      pbc=True)
        calc = elses(command=DUMMY_COMMAND)
        calc.write_input(atoms)

        with open("N2.xyz") as f:
            lines = f.readlines()
        # Check atom 1 position
        parts = lines[2].split()
        self.assertAlmostEqual(float(parts[1]), 1.0, places=5)
        self.assertAlmostEqual(float(parts[2]), 2.0, places=5)
        self.assertAlmostEqual(float(parts[3]), 3.0, places=5)


class TestReadResults(unittest.TestCase):
    """Test read_results parsing of restart.xml and Output.txt."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _write_restart_xml(self, forces_au, natoms):
        """Write a mock restart.xml with forces in atomic units."""
        xml = '<?xml version="1.0"?>\n<root>\n'
        xml += f'  <total_number_of_atoms>{natoms}</total_number_of_atoms>\n'
        for f in forces_au:
            xml += f'  <atom><force>{f[0]} {f[1]} {f[2]}</force></atom>\n'
        xml += '</root>\n'
        with open("restart.xml", "w") as fh:
            fh.write(xml)

    def _write_output_txt(self, energy_ev_per_atom):
        """Write a mock Output.txt matching ELSES output format.
        The parser uses item.split()[4] to extract the energy value,
        so the format is: 'Total Energy summary (eV/atom): <value>'
        where split() yields ['Total', 'Energy', 'summary', '(eV/atom):', '<value>']
        """
        with open("Output.txt", "w") as f:
            f.write("Some header\n")
            f.write(f"  Total Energy summary (eV/atom): {energy_ev_per_atom}\n")
            f.write("Some footer\n")

    def _make_calc(self):
        """Create an elses calculator with results dict initialized."""
        calc = elses(command=DUMMY_COMMAND)
        calc.results = {}
        return calc

    def test_read_results_energy(self):
        """Test energy extraction."""
        natoms = 2
        energy_per_atom = -5.0
        self._write_restart_xml([[0, 0, 0], [0, 0, 0]], natoms)
        self._write_output_txt(energy_per_atom)

        calc = self._make_calc()
        energy, forces = calc.read_results()
        self.assertAlmostEqual(energy, energy_per_atom * natoms, places=5)

    def test_read_results_forces(self):
        """Test force extraction and unit conversion."""
        forces_au = [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0]]
        self._write_restart_xml(forces_au, 2)
        self._write_output_txt(-3.0)

        calc = self._make_calc()
        energy, forces = calc.read_results()

        expected = np.array(forces_au) * eV_angst
        np.testing.assert_array_almost_equal(forces, expected, decimal=5)

    def test_read_results_forces_shape(self):
        """Test forces array shape."""
        forces_au = [[0.1, 0.2, 0.3], [-0.1, -0.2, -0.3], [0.0, 0.0, 0.0]]
        self._write_restart_xml(forces_au, 3)
        self._write_output_txt(-2.0)

        calc = self._make_calc()
        _, forces = calc.read_results()
        self.assertEqual(forces.shape, (3, 3))

    def test_read_results_no_energy(self):
        """Test that missing energy raises PropertyNotPresent."""
        self._write_restart_xml([[0, 0, 0]], 1)
        with open("Output.txt", "w") as f:
            f.write("No energy line here\n")

        calc = self._make_calc()
        with self.assertRaises(PropertyNotPresent):
            calc.read_results()

    def test_read_results_multiple_energy_lines(self):
        """Test that the last energy value is used (for iterative calculations)."""
        self._write_restart_xml([[0, 0, 0]], 1)
        with open("Output.txt", "w") as f:
            f.write("  Total Energy summary (eV/atom): -1.0\n")
            f.write("  Total Energy summary (eV/atom): -2.0\n")
            f.write("  Total Energy summary (eV/atom): -3.0\n")

        calc = self._make_calc()
        energy, _ = calc.read_results()
        # Should use the last value (-3.0), times natoms (1)
        self.assertAlmostEqual(energy, -3.0, places=5)


if __name__ == "__main__":
    unittest.main()
