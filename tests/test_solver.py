"""
Unit tests for AenetSolver integration with ODAT-SE.

Tests cover:
1. Module import and class hierarchy
2. Solver initialization with odatse.Info
3. Template parameter substitution
4. predict.out parsing (_calc_Rfactor)
5. Full evaluate() with mock predict.x
6. Integration with odatse.Runner
"""

import os
import sys
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import odatse
from AenetSolver.aenet import Solver

# Path to mock predict script
MOCK_PREDICT = str(Path(__file__).parent / "mock_predict.sh")

# Save original directory once at module level
_ORIG_DIR = os.getcwd()


def make_test_dir():
    """Create a temporary directory with required test files."""
    tmpdir = tempfile.mkdtemp()

    # Create template.xsf
    with open(os.path.join(tmpdir, "template.xsf"), "w") as f:
        f.write("ATOMS\n")
        f.write("N             0.0000000000        0.0000000000        0.0000000000\n")
        f.write("N             0.0000000000        0.0000000000        value_01\n")

    # Create predict.in
    with open(os.path.join(tmpdir, "predict.in"), "w") as f:
        f.write("TYPES\n1\nN\n\nNETWORKS\n  N N.5t-5t.ann\n\n")
        f.write("FORCES\n\nFILES\n1\nstructure_data.xsf\n")

    # Create dummy .ann potential file
    with open(os.path.join(tmpdir, "N.5t-5t.ann"), "w") as f:
        f.write("dummy ann potential\n")

    # Copy mock predict script
    mock_dest = os.path.join(tmpdir, "mock_predict.sh")
    shutil.copy(MOCK_PREDICT, mock_dest)
    os.chmod(mock_dest, 0o755)

    return tmpdir


def make_info(tmpdir):
    """Create an odatse.Info object for testing."""
    d = {
        "base": {
            "dimension": 1,
            "root_dir": tmpdir,
            "output_dir": "output",
        },
        "algorithm": {
            "name": "pamc",
            "label_list": ["z"],
            "param": {
                "min_list": [0.5],
                "max_list": [2.0],
                "unit_list": [0.015],
            },
        },
        "solver": {
            "name": "aenet",
            "config": {
                "aenet_exec_file": os.path.join(tmpdir, "mock_predict.sh"),
                "aenet_template_file": os.path.join(tmpdir, "template.xsf"),
                "bulk_output_file": os.path.join(tmpdir, "predict.in"),
                "aenet_ann_potential": os.path.join(tmpdir, "N.5t-5t.ann"),
            },
            "param": {
                "string_list": ["value_01"],
            },
        },
    }
    return odatse.Info(d)


class TestImport(unittest.TestCase):
    """Test that AenetSolver can be imported and has correct class hierarchy."""

    def test_import_solver(self):
        from AenetSolver import Solver as S
        self.assertIsNotNone(S)

    def test_inherits_solver_base(self):
        self.assertTrue(issubclass(Solver, odatse.solver.SolverBase))

    def test_solver_name(self):
        self.assertEqual(Solver._name, "aenet")


class TestSolverInit(unittest.TestCase):
    """Test Solver initialization with odatse.Info."""

    def setUp(self):
        os.chdir(_ORIG_DIR)
        self.tmpdir = make_test_dir()
        self.info = make_info(self.tmpdir)

    def tearDown(self):
        os.chdir(_ORIG_DIR)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_init_success(self):
        solver = Solver(self.info)
        self.assertEqual(solver._name, "aenet")
        self.assertEqual(solver.dimension, 1)
        self.assertEqual(solver.string_list, ["value_01"])

    def test_init_sets_paths(self):
        solver = Solver(self.info)
        self.assertTrue(os.access(solver.path_to_solver, os.X_OK))
        self.assertTrue(solver.template_file.exists())
        self.assertTrue(solver.predict_in_file.exists())
        self.assertTrue(solver.ann_potential_file.exists())

    def test_init_missing_template(self):
        os.remove(os.path.join(self.tmpdir, "template.xsf"))
        with self.assertRaises(odatse.exception.InputError):
            Solver(make_info(self.tmpdir))

    def test_init_missing_predict_in(self):
        os.remove(os.path.join(self.tmpdir, "predict.in"))
        with self.assertRaises(odatse.exception.InputError):
            Solver(make_info(self.tmpdir))

    def test_init_missing_ann_potential(self):
        os.remove(os.path.join(self.tmpdir, "N.5t-5t.ann"))
        with self.assertRaises(odatse.exception.InputError):
            Solver(make_info(self.tmpdir))

    def test_init_wrong_dimension(self):
        d = {
            "base": {
                "dimension": 2,
                "root_dir": self.tmpdir,
                "output_dir": "output",
            },
            "algorithm": {"name": "pamc"},
            "solver": {
                "name": "aenet",
                "config": {
                    "aenet_exec_file": os.path.join(self.tmpdir, "mock_predict.sh"),
                    "aenet_template_file": os.path.join(self.tmpdir, "template.xsf"),
                    "bulk_output_file": os.path.join(self.tmpdir, "predict.in"),
                    "aenet_ann_potential": os.path.join(self.tmpdir, "N.5t-5t.ann"),
                },
                "param": {
                    "string_list": ["value_01"],
                },
            },
        }
        with self.assertRaises(odatse.exception.InputError):
            Solver(odatse.Info(d))


class TestTemplateParsing(unittest.TestCase):
    """Test template parameter substitution."""

    def setUp(self):
        os.chdir(_ORIG_DIR)
        self.tmpdir = make_test_dir()
        self.info = make_info(self.tmpdir)
        self.solver = Solver(self.info)

    def tearDown(self):
        os.chdir(_ORIG_DIR)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_format_params_positive(self):
        xs = np.array([1.12345678])
        result = self.solver._format_params(xs)
        self.assertEqual(len(result), 1)
        # value_01 is 8 chars, so formatted value should be truncated to 8 chars
        self.assertEqual(len(result[0]), len("value_01"))
        self.assertIn("1.12345", result[0])

    def test_format_params_negative(self):
        xs = np.array([-0.5])
        result = self.solver._format_params(xs)
        self.assertEqual(len(result), 1)
        self.assertTrue(result[0].startswith("-"))

    def test_write_input_substitution(self):
        proc_dir = Path(self.tmpdir) / "output" / "0"
        os.makedirs(proc_dir, exist_ok=True)
        self.solver.proc_dir = proc_dir
        os.chdir(proc_dir)

        xs = np.array([1.1])
        fitted_x_list = self.solver._format_params(xs)
        folder = self.solver._setup_work_dir(0, 0)
        self.solver._write_input(fitted_x_list, folder)

        output_file = os.path.join(folder, "structure_data.xsf")
        self.assertTrue(os.path.exists(output_file))
        with open(output_file) as f:
            content = f.read()
        # The placeholder should be replaced
        self.assertNotIn("value_01", content)
        self.assertIn("1.1", content)


class TestCalcRfactor(unittest.TestCase):
    """Test predict.out parsing."""

    def setUp(self):
        os.chdir(_ORIG_DIR)
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        os.chdir(_ORIG_DIR)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_parse_predict_out(self):
        predict_out = os.path.join(self.tmpdir, "predict.out")
        with open(predict_out, "w") as f:
            f.write("""
 Number of atoms   :         2
 Number of species :         1
 Total energy      :    -3.36200000 eV
 Energy/atom       :    -1.68100000 eV
""")
        os.chdir(self.tmpdir)

        tmpdir2 = make_test_dir()
        info = make_info(tmpdir2)
        solver = Solver(info)
        result = solver._calc_Rfactor()
        shutil.rmtree(tmpdir2, ignore_errors=True)

        self.assertAlmostEqual(result, -3.36200000 / 2, places=5)

    def test_parse_predict_out_single_atom(self):
        predict_out = os.path.join(self.tmpdir, "predict.out")
        with open(predict_out, "w") as f:
            f.write("""
 Number of atoms   :         1
 Total energy      :    -5.00000000 eV
""")
        os.chdir(self.tmpdir)

        tmpdir2 = make_test_dir()
        info = make_info(tmpdir2)
        solver = Solver(info)
        result = solver._calc_Rfactor()
        shutil.rmtree(tmpdir2, ignore_errors=True)

        self.assertAlmostEqual(result, -5.0, places=5)


class TestEvaluate(unittest.TestCase):
    """Test full evaluate() with mock predict.x."""

    def setUp(self):
        os.chdir(_ORIG_DIR)
        self.tmpdir = make_test_dir()
        self.info = make_info(self.tmpdir)
        self.solver = Solver(self.info)
        # Setup proc_dir like ODAT-SE does
        self.solver.proc_dir = Path(self.tmpdir) / "output" / "0"
        os.makedirs(self.solver.proc_dir, exist_ok=True)

    def tearDown(self):
        os.chdir(_ORIG_DIR)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_evaluate_returns_float(self):
        os.chdir(self.solver.proc_dir)
        result = self.solver.evaluate(np.array([1.1]), args=(0, 0))
        self.assertIsInstance(result, float)

    def test_evaluate_value(self):
        os.chdir(self.solver.proc_dir)
        result = self.solver.evaluate(np.array([1.1]), args=(0, 0))
        # mock outputs: Total energy = -3.362 eV, 2 atoms => -1.681
        self.assertAlmostEqual(result, -3.36200000 / 2, places=5)

    def test_evaluate_creates_work_dir(self):
        os.chdir(self.solver.proc_dir)
        self.solver.evaluate(np.array([1.1]), args=(5, 0))
        work_dir = self.solver.proc_dir / "Log00000005_00000000"
        self.assertTrue(work_dir.exists())

    def test_evaluate_creates_output_files(self):
        os.chdir(self.solver.proc_dir)
        self.solver.evaluate(np.array([1.1]), args=(1, 0))
        work_dir = self.solver.proc_dir / "Log00000001_00000000"
        self.assertTrue((work_dir / "predict.out").exists())
        self.assertTrue((work_dir / "structure_data.xsf").exists())
        self.assertTrue((work_dir / "predict.in").exists())
        self.assertTrue((work_dir / "N.5t-5t.ann").exists())


class TestRunnerIntegration(unittest.TestCase):
    """Test integration with odatse.Runner."""

    def setUp(self):
        os.chdir(_ORIG_DIR)
        self.tmpdir = make_test_dir()
        self.info = make_info(self.tmpdir)
        self.solver = Solver(self.info)
        # Setup proc_dir
        self.solver.proc_dir = Path(self.tmpdir) / "output" / "0"
        os.makedirs(self.solver.proc_dir, exist_ok=True)

    def tearDown(self):
        os.chdir(_ORIG_DIR)
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_runner_creation(self):
        runner = odatse.Runner(self.solver, self.info)
        self.assertIsNotNone(runner)
        self.assertEqual(runner.solver_name, "aenet")

    def test_runner_submit(self):
        runner = odatse.Runner(self.solver, self.info)
        runner.prepare(self.solver.proc_dir)
        os.chdir(self.solver.proc_dir)
        result = runner.submit(np.array([1.1]), (0, 0))
        self.assertIsInstance(result, float)
        self.assertAlmostEqual(result, -3.36200000 / 2, places=5)


if __name__ == "__main__":
    unittest.main()
