Tutorial: Fullerene Molecular Dynamics (ELSES)
================================================

This tutorial explains how to perform an NVT molecular dynamics simulation
of fullerene (C480) using ASE (Atomic Simulation Environment) and
the ELSES electronic structure calculation code, and how to generate
AENET training data from the resulting trajectory.

Sample files are located in ``sample/fullerene_elses_md/``.

Prerequisites
-------------

- ELSES is installed (``elses``, ``elses-xml-generate``)
- ASE is installed
- The ElsesCalculator from this package is installed
- A crystal structure file (``C.cif``)

File Structure
--------------

.. code-block:: text

   sample/fullerene_elses_md/
   ├── do.sh                 # Execution script
   ├── test_elses.py         # MD simulation script
   ├── teach_data_make.py    # Training data generation script
   ├── time_energy_plot.py   # Result visualization script
   └── animation.py          # MD trajectory animation generation

Simulation Settings
-------------------

Key settings in ``test_elses.py``:

ELSES Calculator Settings
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ElsesCalculator import elses

   command = 'elses-xml-generate generate.xml C480.xml ; srun elses config.xml > log.txt'
   calc = elses(command=command)
   supercell.calc = calc

The ELSES execution command is specified in ``command``. ELSES runs in the following two stages:

1. ``elses-xml-generate``: Generation of input XML
2. ``elses``: Execution of electronic structure calculation

Structure Preparation
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ase.io import read
   from ase.build import make_supercell

   atoms = read("C.cif")
   atoms.pbc = [1, 1, 1]
   supercell_matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
   supercell = make_supercell(atoms, supercell_matrix)

A C240 structure is read from a CIF file, and a C480 supercell is created by doubling it in the x direction.

MD Settings
~~~~~~~~~~~

.. code-block:: python

   from ase import units
   from ase.md.nvtberendsen import NVTBerendsen
   from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

   dt = 1.0 * units.fs        # Time step: 1 fs
   temperature_K = 300        # Target temperature: 300 K
   nsteps = 100               # Number of steps: 100

   MaxwellBoltzmannDistribution(supercell, temperature_K=300)

   taut = 0.5 * 10 * units.fs
   dyn = NVTBerendsen(supercell, timestep=dt, temperature=temperature_K,
                      taut=taut, trajectory=traj, logfile='md.log')
   dyn.run(nsteps)

- **Ensemble**: NVT (Berendsen thermostat)
- **Time step**: 1 fs (shorter than the QE version for ELSES stability)
- **Total simulation time**: 100 fs

How to Run
----------

.. code-block:: bash

   cd sample/fullerene_elses_md
   sh do.sh

Output files:

- ``md.log``: Time, energy, kinetic energy, and temperature at each step
- ``benzene_optimization.traj``: ASE trajectory file

Result Visualization
--------------------

.. code-block:: bash

   python3 time_energy_plot.py

.. figure:: ../../common/img/time_energy_elses.png
   :width: 80%
   :align: center

   Time evolution of kinetic energy (eV/atom) for the fullerene supercell calculated by ELSES.
   The red dashed line indicates the theoretical value at 300 K.

Training Data Generation
------------------------

Training data in XSF format for AENET is generated from the MD trajectory.

.. code-block:: bash

   python3 teach_data_make.py

``teach_data_make.py`` performs the following:

1. Reads the structure of each frame from the ``.traj`` file
2. Extracts energy, atomic coordinates, and forces
3. Outputs in XSF format to the ``teach_data/`` directory

.. code-block:: python

   from ase.io import Trajectory

   structures = Trajectory('benzene_optimization.traj')
   for i, atoms in enumerate(structures):
       total_energy = atoms.get_total_energy()
       positions = atoms.get_positions()
       forces = atoms.get_forces()
       # Output in XSF format

Format of the output XSF files:

.. code-block:: text

   # total energy = -XXXX.XX eV

   CRYSTAL
   PRIMVEC
     XX.XXXXXXXX  0.00000000  0.00000000
     0.00000000  XX.XXXXXXXX  0.00000000
     0.00000000  0.00000000  XX.XXXXXXXX
   PRIMCOORD
   480 1
   C  x.xxxxxxxx  y.yyyyyyyy  z.zzzzzzzz  fx.xxxxxxxx  fy.yyyyyyyy  fz.zzzzzzzz
   ...

- For periodic boundary conditions, the ``CRYSTAL`` + ``PRIMVEC`` + ``PRIMCOORD`` format is used
- For isolated systems, the ``ATOMS`` format is used
- Each line: element symbol, x, y, z coordinates, fx, fy, fz forces

The generated training data can be used to build AENET potentials following the procedure in :doc:`tutorial_training`.

Trajectory Animation
--------------------

``animation.py`` can be used to create an animated GIF visualizing the MD trajectory.

.. code-block:: bash

   python3 animation.py --traj benzene_optimization.traj --output animation.gif \
       --natoms 480 --dt 10.0 --nframes 10

The output GIF displays the structural changes and the time evolution of kinetic energy side by side.

.. note::

   ``imagemagick`` is required for animation generation.
   On macOS: ``brew install imagemagick``

Comparison with QE
------------------

Comparison of MD simulations using ELSES and Quantum ESPRESSO:

.. list-table::
   :header-rows: 1
   :widths: 25 35 40

   * - Item
     - QE (:doc:`tutorial_fullerene_qe`)
     - ELSES (this tutorial)
   * - Method
     - DFT (PAW)
     - Tight-binding based
   * - Time step
     - 2 fs
     - 1 fs
   * - Computational cost
     - High (requires 64-core parallelization)
     - Moderate
   * - Accuracy
     - First-principles accuracy
     - Approximate
   * - Use case
     - High-accuracy training data generation
     - Large-scale systems / long-time MD
