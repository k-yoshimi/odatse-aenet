# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

ODAT-SE AENET solver module: integrates AENET machine learning potentials with the ODAT-SE (Open Data Analysis Tool for Science and Engineering) framework for parameter optimization. Also includes an ASE calculator wrapper for ELSES.

ODAT-SE source: `/Users/k-yoshimi/Dropbox/PycharmProjects/ODAT-SE`
Solver-template reference: `odat-se-gallery/data/tutorial/solver-template/`

## Repository Structure

```
src/                          # Source code
  main.py                     # Entry point (python3 src/main.py input.toml)
  AenetSolver/                # ODAT-SE solver for AENET
    __init__.py
    _main.py                  # odatse.initialize() + algorithm dispatch
    aenet.py                  # Solver(odatse.solver.SolverBase)
  ElsesCalculator/            # ASE FileIOCalculator for ELSES
    __init__.py
    elses.py
sample/                       # Tutorial examples
  aenet_pamc/                 # N2 dimer PAMC optimization
  aenet_training/             # AENET training pipeline (generate/train/relax/predict)
  fullerene_qe_md/            # Fullerene MD with Quantum ESPRESSO
  fullerene_elses_md/         # Fullerene MD with ELSES
```

## Install and Run

```bash
pip install .                          # Install package
odatse-aenet input.toml                # CLI
mpiexec -np 4 odatse-aenet input.toml  # With MPI
python3 src/main.py input.toml         # Direct execution
```

## ODAT-SE Solver Pattern

```python
import odatse

class Solver(odatse.solver.SolverBase):
    def __init__(self, info: odatse.Info):
        super().__init__(info)
        self._name = "aenet"

    def evaluate(self, xs, args=(), nprocs=1, nthreads=1) -> float:
        # args = (step, iset); return objective value
```

Main script:
```python
info, run_mode = odatse.initialize()
alg_module = odatse.algorithm.choose_algorithm(info.algorithm["name"])
solver = Solver(info)
runner = odatse.Runner(solver, info)
alg = alg_module.Algorithm(info, runner, run_mode=run_mode)
result = alg.main()
```

## Key Dependencies

- `odatse` (ODAT-SE >= 3.0) -- optimization framework
- `ase` -- atomic simulation environment (for MD samples and ELSES calculator)
- `numpy` -- numerical computation
- External: `predict.x` (AENET), `pw.x` (Quantum ESPRESSO), ELSES

## Configuration

- `input.toml` -- ODAT-SE config with `[base]`, `[solver]`, `[algorithm]` sections
- `template.xsf` -- structure template with placeholders (e.g., `value_01`)
- `predict.in` -- AENET prediction config referencing `.ann` potential files

## Legacy

Old chapter-based structure (`chap3/`, `chap4/`, `chap5/`) can be removed once migration to `src/` + `sample/` is verified.
