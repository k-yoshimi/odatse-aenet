"""
Unit tests for md.log parsing logic used in time_energy_plot.py scripts.

Tests cover:
1. MD log file parsing (time and kinetic energy extraction)
2. Boltzmann constant reference energy calculation
3. Edge cases (empty file, header-only file)
"""

import os
import shutil
import tempfile
import unittest


def parse_md_log(filepath, natoms=480):
    """
    Parse md.log and return time [fs] and ekin [eV/atom] lists.
    This replicates the logic from time_energy_plot.py.
    """
    time = []
    ekin = []
    with open(filepath, "r") as f:
        for line in f:
            if "Time[ps]" in line or line.strip() == "":
                continue
            cols = line.split()
            time.append(float(cols[0]) * 1000)  # ps -> fs
            ekin.append(float(cols[3]) / natoms)  # total Ekin -> per atom
    return time, ekin


class TestMdLogParsing(unittest.TestCase):
    """Test MD log file parsing."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_parse_normal_log(self):
        """Test parsing a typical md.log file."""
        log_content = (
            "Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]\n"
            "0.0000          -100.0       -90.0        10.0       300.0\n"
            "0.0020          -100.1       -89.5        10.6       310.0\n"
            "0.0040          -100.2       -89.8        10.4       305.0\n"
        )
        with open("md.log", "w") as f:
            f.write(log_content)

        time, ekin = parse_md_log("md.log", natoms=2)

        self.assertEqual(len(time), 3)
        self.assertEqual(len(ekin), 3)
        # Time conversion: 0.0000 ps -> 0.0 fs
        self.assertAlmostEqual(time[0], 0.0, places=2)
        # 0.0020 ps -> 2.0 fs
        self.assertAlmostEqual(time[1], 2.0, places=2)
        # Ekin per atom: 10.0 / 2 = 5.0
        self.assertAlmostEqual(ekin[0], 5.0, places=5)

    def test_parse_empty_lines(self):
        """Test that empty lines are skipped."""
        log_content = (
            "Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]\n"
            "\n"
            "0.0000          -100.0       -90.0        10.0       300.0\n"
            "\n"
            "0.0020          -100.1       -89.5        10.6       310.0\n"
        )
        with open("md.log", "w") as f:
            f.write(log_content)

        time, ekin = parse_md_log("md.log", natoms=1)
        self.assertEqual(len(time), 2)

    def test_parse_header_only(self):
        """Test parsing a log with only the header."""
        log_content = "Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]\n"
        with open("md.log", "w") as f:
            f.write(log_content)

        time, ekin = parse_md_log("md.log")
        self.assertEqual(len(time), 0)
        self.assertEqual(len(ekin), 0)

    def test_parse_single_step(self):
        """Test parsing a log with a single timestep."""
        log_content = (
            "Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]\n"
            "0.0010          -200.0       -180.0       20.0       350.0\n"
        )
        with open("md.log", "w") as f:
            f.write(log_content)

        time, ekin = parse_md_log("md.log", natoms=480)
        self.assertEqual(len(time), 1)
        self.assertAlmostEqual(time[0], 1.0, places=2)  # 0.001 ps = 1 fs
        self.assertAlmostEqual(ekin[0], 20.0 / 480, places=8)


class TestBoltzmannReference(unittest.TestCase):
    """Test the kinetic energy reference from Boltzmann distribution."""

    def test_ekin_reference_300K(self):
        """Test Ekin reference at 300 K."""
        k_b = 8.617333e-5  # eV/K
        T = 300
        ekin_ref = (3 / 2) * k_b * T
        # Expected: ~0.03878 eV/atom
        self.assertAlmostEqual(ekin_ref, 0.03878, places=4)

    def test_ekin_reference_scales_with_temperature(self):
        """Test that reference Ekin scales linearly with temperature."""
        k_b = 8.617333e-5
        ekin_300 = (3 / 2) * k_b * 300
        ekin_600 = (3 / 2) * k_b * 600
        self.assertAlmostEqual(ekin_600 / ekin_300, 2.0, places=5)


if __name__ == "__main__":
    unittest.main()
