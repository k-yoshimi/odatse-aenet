Tutorial: AENET Training Pipeline
======================================

In this tutorial, we use the N2 dimer as an example to explain the complete procedure
for creating training data from first-principles calculations (Quantum ESPRESSO)
and building a machine learning potential with AENET.

The sample files are located in ``sample/aenet_training/``.

Workflow Overview
------------------

.. code-block:: text

   1. Structure Generation (relax/)
      └→ Generate random N2 dimer structures
      └→ Perform structural relaxation and energy calculation with QE

   2. Training Data Creation (relax/)
      └→ Extract energies and forces from QE output
      └→ Generate training data in XSF format

   3. Fingerprint Generation (generate/)
      └→ Compute atomic descriptors with AENET generate.x

   4. Neural Network Training (train/)
      └→ Train the potential with AENET train.x

   5. Prediction and Validation (predict/)
      └→ Predict energies of test structures with AENET predict.x
      └→ Compare with QE calculations

Step 1: Structure Generation
-----------------------------

Work in the ``sample/aenet_training/relax/`` directory.

``structure_make.py`` randomly generates N-N bond distances in the range of 0.5 to 2.0 Angstrom
and creates Quantum ESPRESSO input files.

.. code-block:: python

   import random
   structure_list = []
   for i in range(20):
       structure_list.append(random.uniform(0.5, 2.0))

The template file (``template.txt``) is in QE input format, where ``value_01`` is
the placeholder for the N-N bond distance:

.. code-block:: text

   &CONTROL
     calculation = 'relax'
     pseudo_dir  = './'
     outdir      = './'
   /
   &SYSTEM
     ntyp = 1, nat = 2, ibrav = 0
     ecutwfc = 44, ecutrho = 320
     occupations = 'smearing'
     smearing    = 'mp'
     degauss     = 0.01
   /
   ...
   ATOMIC_POSITIONS angstrom
     N 0.00 0.00 0.00
     N 0.00 0.00 value_01

Execution:

.. code-block:: bash

   python3 structure_make.py

This generates ``directory_0/`` through ``directory_19/``, each containing a QE input file.

Run QE in each directory:

.. code-block:: bash

   srun pw.x < n2_dimer.pwi > n2_dimer.pwo

Step 2: Training Data Creation
-------------------------------

Extract energies and interatomic forces from the QE output files and convert them to AENET's XSF format.

``teach_data_make.py`` performs the following:

- Extracts energies from ``.pwo`` files (Ry to eV conversion: ``ase.units.Ry``)
- Extracts atomic coordinates and forces (force conversion: ``ase.units.Ry / 0.529177`` [Ry/Bohr to eV/Angstrom])
- Outputs in XSF format

.. code-block:: bash

   cd directory_0
   python3 ../teach_data_make.py --input n2_dimer.pwo --output-dir teach_data

Output example (teach_data_1.xsf):

.. code-block:: text

   # total energy = -767.51 eV

   ATOMS
   N   0.00000000   0.00000000   0.00000000   0.00000000   0.00000000  -12.34567890
   N   0.00000000   0.00000000   1.04970000   0.00000000   0.00000000   12.34567890

Step 3: Fingerprint Generation
-------------------------------

Work in the ``sample/aenet_training/generate/`` directory.

In ``generate.in``, specify the paths to the training data and the atomic species:

.. code-block:: text

   OUTPUT N2.train

   TYPES
   1
   N 389.83

   SETUPS
   N N.fingerprint.stp

   FILES
   20
   ../relax/directory_0/teach_data/teach_data_1.xsf
   ../relax/directory_1/teach_data/teach_data_1.xsf
   ...

Execution:

.. code-block:: bash

   generate.x generate.in > generate.out

Output:

- ``N2.train``: Binary dataset for training
- ``N.fingerprint.stp``: Fingerprint configuration file

Step 4: Neural Network Training
---------------------------------

Work in the ``sample/aenet_training/train/`` directory.

Set the training parameters in ``train.in``:

.. code-block:: text

   TRAININGSET N2.train
   TESTPERCENT 10
   ITERATIONS  3000

   bfgs

   NETWORKS
   # atom   network           hidden
   # types  file-name         layers   nodes:activation
   N        N.5t-5t.ann       2        10:tanh 10:tanh

Key settings:

- **TESTPERCENT 10**: Use 10% of the training data as test data
- **ITERATIONS 3000**: Run 3000 iterations of BFGS optimization
- **Network architecture**: Input layer -> 10 nodes (tanh) -> 10 nodes (tanh) -> Output layer

Execution:

.. code-block:: bash

   train.x train.in > train.out

Example of training results (from train.out):

.. code-block:: text

   Number of training structures :         18
   Number of testing structures  :          2
   Total number of weights       :        221
   Atomic energy shift           :   -747.530391 eV

Output:

- ``N.5t-5t.ann``: Trained ANN potential file

Step 5: Prediction and Validation
-----------------------------------

Work in the ``sample/aenet_training/predict/`` directory.

Predict energies using the trained potential for test structures
(201 structures with N-N distances from 0.00 to 2.00 Angstrom in 0.01 Angstrom increments).

.. code-block:: bash

   predict.x predict.in > predict.out

Visualization of Results
~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``plot_distance_energy.py`` to create a plot comparing the predicted energies from the ANN potential
with the first-principles calculation results from QE:

.. code-block:: bash

   python3 plot_distance_energy.py \
       --predict-out predict.out \
       --train-out ../train/train.out \
       --qe-dir ../relax \
       --n-structures 20 \
       --output distance_energy_plot.png

This script plots the following three datasets:

- **ANN prediction** (blue): Interatomic distance vs. energy from ``predict.out``
- **QE training data** (red): Results from the relaxation calculation XSF files
- **QE test data** (green): Data points included in the test set from ``train.out``

Calculation Results
~~~~~~~~~~~~~~~~~~~~

The following graph compares the predicted energies from the ANN potential with the first-principles calculation results from QE.

.. figure:: ../../common/img/QE_ANN.png
   :width: 80%
   :align: center

   Energy-distance curve for the N2 dimer. Blue: ANN prediction, Red: QE training data, Green: QE test data.

Key results:

- **ANN most stable distance**: Approximately 1.11 Angstrom (energy: -767.81 eV)
- **QE most stable distance**: Approximately 1.05 Angstrom (energy: -767.51 eV)
- The ANN potential reproduces the QE results well within the range of the training data (20 points)
