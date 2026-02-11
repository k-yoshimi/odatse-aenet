# odatse-aenet

[![Tests](https://github.com/k-yoshimi/odatse-aenet/actions/workflows/test.yml/badge.svg)](https://github.com/k-yoshimi/odatse-aenet/actions/workflows/test.yml)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://k-yoshimi.github.io/odatse-aenet/)

AENET (Atom-centered Neural Network) solver module for [ODAT-SE](https://github.com/issp-center-dev/ODAT-SE) (Open Data Analysis Tool for Science and Engineering).

This package provides:

- **AenetSolver** --- ODAT-SE solver that drives AENET `predict.x` for parameter optimization via grid search, Nelder-Mead, Bayesian optimization, replica exchange, and other algorithms.

## Requirements

- Python >= 3.9
- numpy >= 1.22
- [ODAT-SE](https://github.com/issp-center-dev/ODAT-SE) >= 3.0
- [AENET](http://ann.atomistic.net/) (predict.x, generate.x, train.x)

Optional (for training tutorial):

- [Quantum ESPRESSO](https://www.quantum-espresso.org/) (pw.x)

## Installation

```bash
pip install odat-se[all]       # Install ODAT-SE
pip install .                  # Install this package
```

After installation, the `odatse-aenet` command becomes available:

```bash
odatse-aenet --version
```

See the [documentation](#documentation) for detailed installation instructions including AENET and QE.

## Quick Start

### CLI

```bash
odatse-aenet input.toml
```

### MPI parallel

```bash
mpiexec -np 4 odatse-aenet input.toml
```

## Tutorials

| Directory | Description | Requirements |
|-----------|-------------|--------------|
| `sample/aenet_training/` | AENET training pipeline: structure generation → QE relaxation → fingerprint → training → prediction | QE, AENET |
| `sample/aenet_mapper/` | N2 dimer grid search using AENET + mapper | AENET, ODAT-SE |
| `sample/aenet_minsearch/` | N2 dimer optimization using AENET + Nelder-Mead | AENET, ODAT-SE |

Each tutorial includes a `run_all.sh` script to execute the full workflow.

## Project Structure

```
src/
  AenetSolver/               # ODAT-SE solver for AENET
    aenet.py                 # Solver class (extends odatse.solver.SolverBase)
    _main.py                 # CLI entry point
sample/                      # Tutorial samples
docs/
  ja/                        # Japanese documentation
  en/                        # English documentation
tests/                       # Unit tests (52 tests)
```

## Documentation

Full documentation is available online:

- [**English**](https://k-yoshimi.github.io/odatse-aenet/en/)
- [**日本語**](https://k-yoshimi.github.io/odatse-aenet/ja/)

To build locally:

```bash
cd docs
sh make.sh       # Build both ja/en HTML docs
```

## Developers

- **谷田秀哉** (Shuya Tanida), Tottori University --- Initial code development
- **星健夫** (Takeo Hoshi), National Institute for Fusion Science --- ODAT-SE co-development, research supervision
- **吉見一慶** (Kazuyoshi Yoshimi), ISSP, The University of Tokyo --- Code organization and packaging

## Citation

If you use ODAT-SE in your research, please cite:

> Y. Motoyama, K. Yoshimi, I. Mochizuki, H. Iwamoto, H. Ichinose, and T. Hoshi,
> "Data-analysis software framework 2DMAT and its application to experimental measurements for two-dimensional material structures",
> Computer Physics Communications **280**, 108465 (2022).

## License

[MPL-2.0](https://www.mozilla.org/en-US/MPL/2.0/)
