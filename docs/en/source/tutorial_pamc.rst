Tutorial: Optimization of N2 Dimer Using AENET + PAMC
=====================================================

This tutorial explains how to search for the most stable structure of an N2 dimer
by combining a machine learning potential built with AENET and the
Population Annealing Monte Carlo (PAMC) method in ODAT-SE.

The sample files are located in ``sample/aenet_pamc/``.

Prerequisites
-------------

- AENET's ``predict.x`` must be installed
- A trained ANN potential file (``N.5t-5t.ann``) must be prepared
  (see :doc:`tutorial_training`)

File Structure
--------------

.. code-block:: text

   sample/aenet_pamc/
   ├── do.sh                  # Execution script
   ├── input.toml             # ODAT-SE configuration file
   ├── predict.in             # AENET prediction settings
   ├── template.xsf           # Structure template
   ├── plot_histogram.py      # Histogram plot of results
   ├── summarize_results.py   # Summary of results per temperature
   └── batch_plot.py          # Batch histogram plotting

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
   name = "pamc"
   label_list = ["z"]

   [algorithm.param]
   min_list = [0.5]
   max_list = [2.0]
   unit_list = [0.015]

   [algorithm.pamc]
   numsteps_per_T = 40
   downTsteps_resampling = 2
   Tnum = 41
   Tmin = 1
   Tmax = 1000
   replica_per_prop = 10
   Tlogspace = true
   resampling_bool = true
   use_so_for_resampling = true

PAMC algorithm settings:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``min_list`` / ``max_list``
     - Search range: N-N distance 0.5 to 2.0 Angstrom
   * - ``unit_list``
     - Unit step size for Monte Carlo steps: 0.015 Angstrom
   * - ``Tnum``
     - Number of temperature points: 41
   * - ``Tmin`` / ``Tmax``
     - Temperature range: 1 to 1000 K
   * - ``Tlogspace``
     - Divide temperatures on a logarithmic scale
   * - ``numsteps_per_T``
     - Number of Monte Carlo steps per temperature: 40
   * - ``replica_per_prop``
     - Number of replicas per proposal: 10
   * - ``resampling_bool``
     - Enable resampling

Template File
-------------

``template.xsf`` is the structure template for the N2 dimer:

.. code-block:: text

   ATOMS
   N             0.0000000000        0.0000000000        0.0000000000
   N             0.0000000000        0.0000000000        value_01

``value_01`` is the N-N bond distance, which PAMC explores in the range of 0.5 to 2.0 Angstrom.

How to Run
----------

.. code-block:: bash

   cd sample/aenet_pamc
   sh do.sh

Or run directly:

.. code-block:: bash

   odatse-aenet input.toml

MPI parallel execution:

.. code-block:: bash

   mpiexec -np 2 odatse-aenet input.toml

Output
------

The calculation results are generated in the ``output/`` directory:

- ``output/0/`` , ``output/1/`` , ...: Output for each MPI rank
- Each rank contains files such as ``ColorMap.txt``, ``BestSteps.txt``, etc.
- Working directories ``Log{step}_{iset}/`` contain the solver input/output files

Results
-------

As a result of PAMC, histograms of the N-N bond distance at each temperature are obtained.

At low temperatures (T -> 1 K), the distribution sharply concentrates around the energy minimum (approximately 1.1 Angstrom).
At high temperatures (T -> 1000 K), the distribution becomes broad and flat over a wide range.

.. note::

   After execution, you can visualize the probability distribution of the N-N bond distance
   at each temperature by creating histograms from the output data.

Analysis of Results
~~~~~~~~~~~~~~~~~~~

To check the optimal parameters from the output files:

.. code-block:: bash

   # Find the step with lowest R-factor
   sort -k2 -n output/0/ColorMap.txt | head -5

This yields the N-N bond distance at the energy minimum (approximately 1.1 Angstrom) as the optimal solution.
This value is in good agreement with the experimental bond distance of the N2 molecule (1.098 Angstrom).

Analysis Tools
~~~~~~~~~~~~~~

Scripts are provided for analyzing and visualizing the PAMC output data.

**1. Summarizing Results (summarize_results.py)**

Aggregates weight data from ``weight.txt`` for each temperature step
and outputs normalized/unnormalized weight files:

.. code-block:: bash

   python3 summarize_results.py -i input.toml -d my_data --idnum 0

**2. Histogram Plot (plot_histogram.py)**

Creates a weighted 1D histogram for a specified variable from the summary file:

.. code-block:: bash

   python3 plot_histogram.py <summary-file> <temperature-step> tau x1

**3. Batch Plotting (batch_plot.py)**

Generates histograms for all temperature steps at once:

.. code-block:: bash

   python3 batch_plot.py plot_histogram.py <summary-directory> --variable x1
