Tutorial: Minimization with AENET + minsearch
=====================================================

This tutorial explains how to find the optimal N-N bond distance of an N2 dimer
by combining a machine learning potential built with AENET and the
minsearch (Nelder-Mead method) algorithm in ODAT-SE.

The sample files are located in ``sample/aenet_minsearch/``.

Prerequisites
-------------

- AENET's ``predict.x`` must be installed
- A trained ANN potential file (``N.5t-5t.ann``) must be prepared
  (see :doc:`tutorial_training`)

File Structure
--------------

.. code-block:: text

   sample/aenet_minsearch/
   ├── input.toml             # ODAT-SE configuration file
   ├── predict.in             # AENET prediction settings
   ├── template.xsf           # Structure template
   └── run_all.sh             # Execution script

Description of input.toml
--------------------------

.. code-block:: toml

   [base]
   dimension = 1
   output_dir = "output"

Sets the dimension of the parameter space to 1 (N-N bond distance only).

.. code-block:: toml

   [solver]
   name = "aenet"

   [solver.config]
   aenet_exec_file = "predict.x"
   aenet_ann_potential = "N.5t-5t.ann"

   [solver.param]
   string_list = ["value_01"]

Specifies AENET as the solver, along with the path to ``predict.x`` and the ANN potential file.
``string_list`` defines the placeholders in the template.

.. code-block:: toml

   [algorithm]
   name = "minsearch"
   label_list = ["z"]

   [algorithm.param]
   min_list = [0.5]
   max_list = [2.0]
   initial_list = [1.25]

   [algorithm.minimize]
   initial_scale_list = [0.25]
   xatol = 0.0001
   fatol = 0.0001
   maxiter = 100

Minsearch algorithm settings:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``min_list`` / ``max_list``
     - Search range: N-N distance 0.5 to 2.0 Angstrom
   * - ``initial_list``
     - Initial value: 1.25 Angstrom
   * - ``initial_scale_list``
     - Scale of the initial simplex: 0.25 Angstrom
   * - ``xatol``
     - Convergence threshold for parameters: 0.0001 Angstrom
   * - ``fatol``
     - Convergence threshold for function value: 0.0001 eV/atom
   * - ``maxiter``
     - Maximum number of iterations: 100

Template File
-------------

``template.xsf`` is the structure template for the N2 dimer:

.. code-block:: text

   ATOMS
   N             0.0000000000        0.0000000000        0.0000000000
   N             0.0000000000        0.0000000000        value_01

``value_01`` is the N-N bond distance to be optimized by minsearch.

How to Run
----------

Using the execution script:

.. code-block:: bash

   cd sample/aenet_minsearch
   sh run_all.sh

Or run directly:

.. code-block:: bash

   odatse-aenet input.toml

The Nelder-Mead optimization converges in 13 iterations (26 function evaluations) and completes in about 1 second.

Output
------

The calculation results are generated in the ``output/`` directory:

- ``output/res.txt``: Optimization result (optimal parameters and objective function value)
- ``output/SimplexData.txt``: Iteration history of the Nelder-Mead simplex

``res.txt`` contains the optimized N-N bond distance and its corresponding energy (eV/atom).

Format of ``SimplexData.txt``:

.. code-block:: text

   # step  z  fx

Each line represents the parameter values and function value at each step of the simplex.

Results
-------

The Nelder-Mead optimization yields an optimal N-N bond distance of approximately 1.1 Angstrom.

.. figure:: ../../common/img/minsearch_convergence.png
   :width: 500px
   :align: center

   Convergence of the N-N bond distance by the Nelder-Mead method. Top: N-N distance. Bottom: Energy (eV/atom).

To check the optimization result:

.. code-block:: bash

   cat output/res.txt
