AENET Solver
============

Overview
--------

AenetSolver is a solver class that inherits from ``odatse.solver.SolverBase`` in ODAT-SE.
It runs AENET's ``predict.x`` as a subprocess, substitutes parameters into a structure template,
and calculates interatomic energies.

How the Solver Works
--------------------

When the ``evaluate`` method is called, the following steps are performed:

1. Create a working directory ``Log{step}_{iset}/``
2. Copy ``predict.in`` and the ANN potential file
3. Replace placeholders in the template file (e.g., ``value_01``) with parameter values to generate the input structure file
4. Execute ``predict.x``
5. Extract the total energy and the number of atoms from ``predict.out``
6. Return the energy per atom (eV/atom) as the R-factor

Source Code Structure
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

   src/AenetSolver/
     __init__.py        # Export of the Solver class
     _main.py           # CLI entry point (odatse.initialize + algorithm dispatch)
     aenet.py           # Solver class implementation

Input Files
-----------

input.toml
~~~~~~~~~~

This is the standard ODAT-SE TOML configuration file.

.. code-block:: toml

   [base]
   dimension = 1
   output_dir = "output"

   [solver]
   name = "aenet"

   [solver.config]
   aenet_exec_file = "predict.x"           # Path to predict.x
   aenet_ann_potential = "N.5t-5t.ann"      # ANN potential file
   aenet_template_file = "template.xsf"    # Template structure file (default)
   aenet_input_file = "structure_data.xsf"  # Generated input file name (default)
   bulk_output_file = "predict.in"          # AENET configuration file (default)

   [solver.param]
   string_list = ["value_01"]              # Placeholders in the template

   [solver.post]
   remove_work_dir = false                 # Whether to remove the working directory after calculation

   [algorithm]
   name = "mapper"                          # Algorithm to use
   label_list = ["z"]

   [algorithm.param]
   min_list = [0.8]
   max_list = [1.4]
   num_list = [101]

solver.config Section
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``aenet_exec_file``
     - ``predict.x``
     - Path to the AENET predict.x executable
   * - ``aenet_ann_potential``
     - ``N.10t-10t.ann``
     - Pre-trained ANN potential file
   * - ``aenet_template_file``
     - ``template.xsf``
     - Structure template file
   * - ``aenet_input_file``
     - ``structure_data.xsf``
     - Name of the generated structure input file
   * - ``bulk_output_file``
     - ``predict.in``
     - AENET configuration file

solver.param Section
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - Parameter
     - Default
     - Description
   * - ``string_list``
     - ``["value_01", "value_02"]``
     - List of placeholders to be replaced in the template

Template File
~~~~~~~~~~~~~

The template file (``template.xsf``) is a structure file in XSF format,
where placeholders are placed at the positions of the optimization parameters.

.. code-block:: text

   ATOMS
   N             0.0000000000        0.0000000000        0.0000000000
   N             0.0000000000        0.0000000000        value_01

In this example, ``value_01`` corresponds to the N-N bond distance and is optimized by the algorithm.

predict.in
~~~~~~~~~~

This is the configuration file for AENET prediction calculations.

.. code-block:: text

   TYPES
   1
   N

   NETWORKS
     N N.5t-5t.ann

   FORCES

   FILES
   1
   structure_data.xsf

How to Run
----------

CLI execution:

.. code-block:: bash

   odatse-aenet input.toml

MPI parallel execution:

.. code-block:: bash

   mpiexec -np 4 odatse-aenet input.toml

Direct execution from a script:

.. code-block:: bash

   python3 src/main.py input.toml

