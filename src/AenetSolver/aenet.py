from typing import List, Tuple
import itertools
import os
import shutil
import subprocess
import re
from pathlib import Path

import numpy as np

import odatse
from odatse import exception


class Solver(odatse.solver.SolverBase):
    """AENET solver for ODAT-SE.

    Wraps the AENET ``predict.x`` executable as an ODAT-SE solver.
    Substitutes parameter values into a template XSF structure file,
    runs the AENET prediction, and returns the total energy per atom.

    Attributes
    ----------
    path_to_solver : Path
        Path to the ``predict.x`` executable.
    string_list : list of str
        Placeholder strings in the template file to be replaced.
    template_file : Path
        Path to the template XSF structure file.
    predict_in_file : Path
        Path to the ``predict.in`` configuration file.
    ann_potential_file : Path
        Path to the ANN potential file.
    input_filename : str
        Name of the output structure file (default: ``structure_data.xsf``).
    remove_work_dir : bool
        Whether to remove working directories after evaluation.
    """

    path_to_solver: Path
    _name = "aenet"

    def __init__(self, info: odatse.Info):
        """Initialize the AENET solver.

        Parameters
        ----------
        info : odatse.Info
            ODAT-SE configuration object containing ``solver.config`` and
            ``solver.param`` sections.

        Raises
        ------
        odatse.exception.InputError
            If required files are missing or ``len(string_list) != dimension``.
        """
        super().__init__(info)
        self._name = "aenet"

        info_config = info.solver.get("config", {})

        # Locate predict.x executable
        p2solver = info_config.get("aenet_exec_file", "predict.x")
        if os.path.dirname(p2solver) != "":
            self.path_to_solver = self.root_dir / Path(p2solver).expanduser()
        else:
            for search_dir in itertools.chain([self.root_dir], os.environ["PATH"].split(":")):
                self.path_to_solver = Path(search_dir) / p2solver
                if os.access(self.path_to_solver, mode=os.X_OK):
                    break
        if not os.access(self.path_to_solver, mode=os.X_OK):
            raise exception.InputError(f"ERROR: solver ({p2solver}) is not found")

        # Read solver parameters
        info_param = info.solver.get("param", {})
        v = info_param.setdefault("string_list", ["value_01", "value_02"])
        if len(v) != self.dimension:
            raise exception.InputError(
                f"ERROR: len(string_list) != dimension ({len(v)} != {self.dimension})"
            )
        self.string_list = v

        # Template file for structure input
        self.input_filename = info_config.get("aenet_input_file", "structure_data.xsf")

        self.template_file = self._resolve_path(
            info_config.get("aenet_template_file", "template.xsf"),
            "aenet_template_file",
        )
        self._check_template()

        # predict.in file
        self.predict_in_file = self._resolve_path(
            info_config.get("bulk_output_file", "predict.in"),
            "bulk_output_file",
        )

        # ANN potential file
        self.ann_potential_file = self._resolve_path(
            info_config.get("aenet_ann_potential", "N.10t-10t.ann"),
            "aenet_ann_potential",
        )

        # Post-processing options
        info_post = info.solver.get("post", {})
        self.remove_work_dir = info_post.get("remove_work_dir", False)

    def _resolve_path(self, filename: str, label: str) -> Path:
        """Resolve a config file path relative to root_dir and verify existence.

        Parameters
        ----------
        filename : str
            Filename or relative path from the configuration.
        label : str
            Human-readable label for error messages.

        Returns
        -------
        Path
            Absolute resolved path.

        Raises
        ------
        odatse.exception.InputError
            If the resolved path does not exist.
        """
        path = self.root_dir / Path(filename).expanduser()
        if not path.exists():
            raise exception.InputError(f"ERROR: {label} ({path}) does not exist")
        return path

    def evaluate(self, xs: np.ndarray, args: Tuple = (), nprocs: int = 1, nthreads: int = 1) -> float:
        """
        Evaluate AENET prediction for given parameters.

        Parameters
        ----------
        xs : np.ndarray
            Parameter values (atomic coordinates to substitute).
        args : Tuple, optional
            (step, iset) for work directory naming.
        nprocs : int, optional
            Number of processes.
        nthreads : int, optional
            Number of threads.

        Returns
        -------
        float
            R-factor (total energy per atom).
        """
        step = args[0] if len(args) > 0 else 0
        iset = args[1] if len(args) > 1 else 0

        fitted_x_list = self._format_params(xs)
        folder_name = self._setup_work_dir(step, iset)
        self._write_input(fitted_x_list, folder_name)
        self.work_dir = self.proc_dir / folder_name

        cwd = os.getcwd()
        os.chdir(self.work_dir)
        try:
            command = f"{self.path_to_solver} predict.in > predict.out"
            subprocess.run(command, shell=True, check=True)
            result = self._calc_Rfactor()
        finally:
            os.chdir(cwd)

        if self.remove_work_dir:
            def rmtree_error_handler(function, path, excinfo):
                print(f"WARNING: Failed to remove a working directory, {path}")
            shutil.rmtree(self.work_dir, onerror=rmtree_error_handler)

        return result

    def _format_params(self, xs: np.ndarray) -> List[str]:
        """Format parameter values as fixed-width strings for template substitution.

        Each value is formatted to ``".8f"`` precision and truncated to the
        length of the corresponding placeholder string.

        Parameters
        ----------
        xs : numpy.ndarray
            Parameter values, shape ``(dimension,)``.

        Returns
        -------
        list of str
            Formatted strings, one per parameter.
        """
        fitted_x_list = []
        for index in range(self.dimension):
            fitted_value = " " if xs[index] >= 0 else ""
            fitted_value += format(xs[index], ".8f")
            fitted_value = fitted_value[: len(self.string_list[index])]
            fitted_x_list.append(fitted_value)
        return fitted_x_list

    def _setup_work_dir(self, step: int, iset: int) -> str:
        """Create a work directory and copy required input files.

        Parameters
        ----------
        step : int
            Step index for directory naming.
        iset : int
            Set index for directory naming.

        Returns
        -------
        str
            Name of the created directory (e.g., ``Log00000005_00000000``).
        """
        folder_name = f"Log{step:08d}_{iset:08d}"
        os.makedirs(folder_name, exist_ok=True)
        shutil.copy(self.predict_in_file, os.path.join(folder_name, self.predict_in_file.name))
        shutil.copy(self.ann_potential_file, os.path.join(folder_name, self.ann_potential_file.name))
        return folder_name

    def _write_input(self, fitted_x_list: List[str], folder_name: str) -> None:
        """Write the structure input file with parameter substitution.

        Reads the template file, replaces each placeholder string with
        the corresponding formatted value, and writes the result into
        the work directory.

        Parameters
        ----------
        fitted_x_list : list of str
            Formatted parameter values from :meth:`_format_params`.
        folder_name : str
            Name of the work directory.
        """
        output_path = os.path.join(folder_name, self.input_filename)
        with open(self.template_file, "r") as fin, open(output_path, "w") as fout:
            for line in fin:
                for index in range(self.dimension):
                    if self.string_list[index] in line:
                        line = line.replace(self.string_list[index], fitted_x_list[index])
                fout.write(line)

    def _check_template(self) -> None:
        """Verify that all placeholder strings exist in the template file.

        Raises
        ------
        odatse.exception.InputError
            If any placeholder from ``string_list`` is not found in the template.
        """
        found = [False] * self.dimension
        with open(self.template_file, "r") as f:
            for line in f:
                for index, placeholder in enumerate(self.string_list):
                    if placeholder in line:
                        found[index] = True
        if not all(found):
            missing = [label for label, f in zip(self.string_list, found) if not f]
            raise exception.InputError(
                "ERROR: the following labels do not appear in the template file:\n"
                + "\n".join(missing)
            )

    def _calc_Rfactor(self) -> float:
        """Parse ``predict.out`` and compute the total energy per atom.

        Returns
        -------
        float
            Total energy divided by the number of atoms (eV/atom).

        Raises
        ------
        RuntimeError
            If total energy or number of atoms cannot be parsed.
        """
        with open("predict.out", "r") as f:
            text = f.read()

        match = re.search(r'Total energy\s+:\s+([-\d.]+) eV', text)
        if match is None:
            raise RuntimeError("Failed to parse total energy from predict.out")
        total_energy = float(match.group(1))

        match = re.search(r'Number of atoms\s+:\s+(\d+)', text)
        if match is None:
            raise RuntimeError("Failed to parse number of atoms from predict.out")
        num_atoms = int(match.group(1))

        return total_energy / num_atoms
