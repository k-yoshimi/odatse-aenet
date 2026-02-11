Tutorial: Ni Bulk Molecular Dynamics with EMT
===============================================

This tutorial explains how to perform an NVT molecular dynamics simulation of
Ni bulk using the EMT (Effective Medium Theory) potential built into ASE.

Since it does not require any external calculation engine (such as QE or ELSES)
and runs entirely within ASE, it is well suited for verifying the basic workflow
of an MD simulation.

The sample files are located in ``sample/emt_md/``.

Prerequisites
-------------

- ASE is installed

File Structure
--------------

.. code-block:: text

   sample/emt_md/
   └── run_md.py   # MD simulation script

Simulation Settings
-------------------

``run_md.py`` performs NVT MD with two-stage temperature control:

1. **Equilibration**: 200 steps at low temperature (10 K)
2. **Production**: 400 steps after raising the temperature to 500 K

Structure Preparation
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from ase.build import bulk, make_supercell

   atoms = bulk("Ni", "fcc", a=3.5, cubic=True)
   atoms = make_supercell(atoms, np.diag([3, 3, 3]))
   atoms.calc = EMT()

The FCC Ni unit cell is expanded into a 3x3x3 supercell (108 atoms) and the
EMT potential is assigned.

MD Settings
~~~~~~~~~~~

.. code-block:: python

   from ase import units
   from ase.md.nvtberendsen import NVTBerendsen
   from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

   dt = 2.0 * units.fs
   MaxwellBoltzmannDistribution(atoms, temperature_K=10)

   dyn = NVTBerendsen(atoms, dt, temperature=10, taut=20*units.fs,
                      trajectory="md.traj")
   dyn.run(200)       # Equilibration

   dyn.set_temperature(500)
   dyn.run(400)       # Production

- **Ensemble**: NVT (Berendsen thermostat)
- **Time step**: 2 fs
- **Thermostat time constant**: 20 fs
- **Initial temperature**: 10 K, raised to 500 K

How to Run
----------

.. code-block:: bash

   cd sample/emt_md
   python3 run_md.py

Each parameter can be changed via command-line options:

.. code-block:: bash

   python3 run_md.py --element Ni --supercell 3 \
       --temp-equil 10 --nsteps-equil 200 \
       --temp-prod 500 --nsteps-prod 400

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Option
     - Default
     - Description
   * - ``--element``
     - ``Ni``
     - Element symbol
   * - ``--supercell``
     - ``3``
     - Supercell multiplier (same in each direction)
   * - ``--temp-equil``
     - ``10``
     - Equilibration temperature [K]
   * - ``--nsteps-equil``
     - ``200``
     - Number of equilibration steps
   * - ``--temp-prod``
     - ``500``
     - Production temperature [K]
   * - ``--nsteps-prod``
     - ``400``
     - Number of production steps
   * - ``--traj``
     - ``md.traj``
     - Output trajectory file

Output
------

- The time and temperature are printed to the console at each step
- ``md.traj``: ASE trajectory file

.. code-block:: text

   --- Equilibration at 10 K ---
   time=    0 fs T= 10 K
   time=   40 fs T= 12 K
   ...
   --- Production at 500 K ---
   time=  400 fs T= 489 K
   time=  440 fs T= 503 K
   ...

After stabilizing near 10 K during the equilibration phase, the temperature is
raised to 500 K, and you can observe the Berendsen thermostat relaxing the
system toward the target temperature.

.. note::

   EMT supports Cu, Ag, Au, Ni, Pd, Pt, C, N, and O.
   To use other elements, switch to an appropriate calculation engine (e.g., QE).
