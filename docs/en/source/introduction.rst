Introduction
============

What is odatse-aenet
---------------------

odatse-aenet is a solver module for integrating the machine learning potential `AENET <http://ann.atomistic.net/>`_ (Atom-centered Neural Network) with
the `ODAT-SE <https://github.com/issp-center-dev/ODAT-SE>`_ (Open Data Analysis Tool for Science and Engineering) framework.

Using the optimization algorithms provided by ODAT-SE, you can perform parameter searches with machine learning potentials constructed by AENET.

Available algorithms
~~~~~~~~~~~~~~~~~~~~~~

The following algorithms are available in the ODAT-SE framework:

- Nelder-Mead method (``minsearch``)
- Grid search method (``mapper``)
- Bayesian optimization (``bayes``)
- Replica exchange Monte Carlo method (``exchange``)
- Population annealing Monte Carlo method (``pamc``)

Package structure
~~~~~~~~~~~~~~~~~~

- **AenetSolver**: Executes AENET predict.x as an ODAT-SE solver class to perform energy calculations.

Developers
----------

- **谷田秀哉** (Graduate School of Sustainability Science, Tottori University) --- Initial code development (Master's thesis)
- **星健夫** (National Institute for Fusion Science) --- ODAT-SE co-development, research supervision
- **吉見一慶** (Institute for Solid State Physics, The University of Tokyo) --- Code organization and packaging

License
-------

This software is released under the Mozilla Public License 2.0 (MPL-2.0).

Citation
--------

When publishing research results using ODAT-SE, please cite the following reference:

  Y. Motoyama, K. Yoshimi, I. Mochizuki, H. Iwamoto, H. Ichinose, and T. Hoshi,
  "Data-analysis software framework 2DMAT and its application to experimental measurements for two-dimensional material structures",
  Computer Physics Communications **280**, 108465 (2022).
