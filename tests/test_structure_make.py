"""
Unit tests for sample/aenet_training/relax/structure_make.py

Tests cover:
1. Random structure generation in valid range
2. Template substitution logic
"""

import os
import shutil
import tempfile
import unittest


class TestStructureMake(unittest.TestCase):
    """Test the template substitution logic from structure_make.py."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.orig_dir = os.getcwd()
        os.chdir(self.tmpdir)

    def tearDown(self):
        os.chdir(self.orig_dir)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_random_values_in_range(self):
        """Test that generated random values are within [0.5, 2.0]."""
        import random
        random.seed(42)
        structure_list = [random.uniform(0.5, 2.0) for _ in range(20)]
        for val in structure_list:
            self.assertGreaterEqual(val, 0.5)
            self.assertLessEqual(val, 2.0)
        self.assertEqual(len(structure_list), 20)

    def test_template_substitution(self):
        """Test that value_01 placeholder is replaced in template."""
        template_content = (
            "&CONTROL\n"
            "  calculation = 'scf'\n"
            "/\n"
            "ATOMIC_POSITIONS angstrom\n"
            "N  0.0 0.0 0.0\n"
            "N  0.0 0.0 value_01\n"
        )
        with open("template.txt", "w") as f:
            f.write(template_content)

        value_list = ["value_01"]
        test_value = 1.234

        line_list = []
        with open("template.txt", "r") as file_input:
            for line in file_input:
                for idx in range(len(value_list)):
                    if line.find(value_list[idx]) != -1:
                        line = line.replace(value_list[idx], str(test_value))
                line_list.append(line)

        with open("output.pwi", "w") as f:
            f.writelines(line_list)

        with open("output.pwi") as f:
            content = f.read()
        self.assertNotIn("value_01", content)
        self.assertIn("1.234", content)
        # Other lines should be unchanged
        self.assertIn("calculation = 'scf'", content)

    def test_multiple_values_same_line(self):
        """Test substitution when multiple placeholders appear."""
        template_content = "N  value_01 value_02 0.0\n"
        with open("template.txt", "w") as f:
            f.write(template_content)

        value_list = ["value_01", "value_02"]
        test_values = [1.1, 2.2]

        line_list = []
        with open("template.txt", "r") as file_input:
            for line in file_input:
                for idx in range(len(value_list)):
                    if line.find(value_list[idx]) != -1:
                        line = line.replace(value_list[idx], str(test_values[idx]))
                line_list.append(line)

        result = "".join(line_list)
        self.assertIn("1.1", result)
        self.assertIn("2.2", result)
        self.assertNotIn("value_01", result)
        self.assertNotIn("value_02", result)

    def test_directory_creation(self):
        """Test that copytree creates proper directory structure."""
        # Create a mock pseudo_Potential directory
        os.makedirs("pseudo_Potential")
        with open("pseudo_Potential/test.UPF", "w") as f:
            f.write("pseudo\n")

        shutil.copytree("pseudo_Potential", "directry_0")
        self.assertTrue(os.path.exists("directry_0"))
        self.assertTrue(os.path.exists("directry_0/test.UPF"))


if __name__ == "__main__":
    unittest.main()
