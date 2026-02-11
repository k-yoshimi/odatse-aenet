Installation
============

Requirements
------------

- Python >= 3.9
- numpy >= 1.22
- `ODAT-SE <https://github.com/issp-center-dev/ODAT-SE>`_ >= 3.0
- `AENET <http://ann.atomistic.net/>`_ (predict.x, generate.x, train.x)

The following are also required depending on the sample:

- `ASE <https://wiki.fysik.dtu.dk/ase/>`_ (Atomic Simulation Environment)
- `Quantum ESPRESSO <https://www.quantum-espresso.org/>`_ (pw.x)
- `ELSES <https://www.elses.jp/>`_

Installing ODAT-SE
-----------------------

To install from PyPI:

.. code-block:: bash

   pip3 install odat-se[all]

To install from source:

.. code-block:: bash

   git clone https://github.com/issp-center-dev/ODAT-SE.git
   cd ODAT-SE
   pip3 install .[all]

Installing AENET
---------------------

AENET (Atomic Energy Network) is used to build machine learning potentials.
Building requires a Fortran compiler (gfortran) and LAPACK/BLAS.

Prerequisites for macOS
~~~~~~~~~~~~~~~~~~~~~~~

If gfortran is not installed, install it via Homebrew:

.. code-block:: bash

   brew install gcc

This makes the ``gfortran`` command available.
On macOS, the Accelerate framework is available by default as a LAPACK/BLAS implementation.

Obtaining the source code
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone https://github.com/atomisticnet/aenet.git
   cd aenet

Building the L-BFGS-B library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before building AENET itself, you need to build the bundled L-BFGS-B library first:

.. code-block:: bash

   cd lib
   make static
   cd ..

Verify that ``lib/liblbfgsb.a`` has been generated.

Building AENET
~~~~~~~~~~~~~~~~~~

In the ``src/`` directory, build by specifying the Makefile appropriate for your environment:

.. code-block:: bash

   cd src
   make -f makefiles/Makefile.gfortran_serial_MacOS   # macOS (Apple Silicon / Intel)
   cd ..

.. note::

   For Linux environments, choose one of the following depending on whether MPI is available:

   - ``Makefile.gfortran_serial`` --- Serial version
   - ``Makefile.gfortran_mpi`` --- MPI parallel version
   - ``Makefile.gfortran_openblas_serial`` --- OpenBLAS version

   See ``src/makefiles/`` for a full list of available Makefiles.

Once the build is complete, the following executables are generated under ``bin/``:

- ``generate.x`` --- Generation of training data (structural descriptors)
- ``train.x`` --- Training of neural network potentials
- ``predict.x`` --- Energy prediction using trained potentials

.. note::

   The executable filenames include version and compiler suffixes
   (e.g., ``predict.x-2.0.4-gfortran_serial``).
   Create symbolic links as needed:

   .. code-block:: bash

      cd bin
      ln -s predict.x-2.0.4-gfortran_serial predict.x
      ln -s generate.x-2.0.4-gfortran_serial generate.x
      ln -s train.x-2.0.4-gfortran_serial train.x

Setting up the PATH:

.. code-block:: bash

   export PATH=/path/to/aenet/bin:$PATH

Verification:

.. code-block:: bash

   generate.x
   # If "generate.x - training set generation" is displayed, the installation was successful

.. note::

   For more details, see the `AENET official website <http://ann.atomistic.net/>`_.

Installing ASE
-------------------

ASE (Atomic Simulation Environment) is used for molecular dynamics simulations
and interfacing with various computational engines.

.. code-block:: bash

   pip3 install ase

Verifying the installation:

.. code-block:: bash

   python3 -c "import ase; print(ase.__version__)"

.. note::

   For more details, see the `ASE official documentation <https://wiki.fysik.dtu.dk/ase/>`_.

Installing Quantum ESPRESSO
---------------------------------

Quantum ESPRESSO is a first-principles calculation engine used in the fullerene MD tutorial.

Prerequisites (macOS)
~~~~~~~~~~~~~~~~~~~~~

cmake is required:

.. code-block:: bash

   brew install cmake

Building from source
~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone --depth 1 --branch qe-7.4 https://github.com/QEF/q-e.git qe-7.4
   cd qe-7.4

.. warning::

   Versions v7.3.1 and earlier fail to build the mbd library with GCC 15 (gfortran 15).
   **Using v7.4 or later is recommended.**

Configure the build with cmake. For the case without MPI and using the internal FFTW:

.. code-block:: bash

   cmake -B build \
     -DCMAKE_Fortran_COMPILER=gfortran \
     -DCMAKE_C_COMPILER=gcc-15 \
     -DQE_ENABLE_MPI=OFF \
     -DQE_FFTW_VENDOR=Internal

.. note::

   - Adjust ``gcc-15`` according to your environment
     (e.g., on Linux where ``gcc`` is GNU, simply use ``gcc``).
   - If you need the MPI parallel version, remove ``-DQE_ENABLE_MPI=OFF``
     and specify ``-DCMAKE_Fortran_COMPILER=mpif90``.
   - If an external FFTW3 is already installed, ``-DQE_FFTW_VENDOR=Internal`` is not needed.

To build only pw.x:

.. code-block:: bash

   cmake --build build --target pw -j4

To build everything:

.. code-block:: bash

   cmake --build build -j4

After building, add ``pw.x`` to your PATH:

.. code-block:: bash

   export PATH=/path/to/qe-7.4/build/bin:$PATH

Verification:

.. code-block:: bash

   pw.x --version
   # If "Program PWSCF v.7.4" is displayed, the installation was successful

.. note::

   - Pseudopotential files (``.UPF``) are required.
     They can be obtained from the
     `SSSP library <https://www.materialscloud.org/discover/sssp/>`_
     and other sources.
   - For more details, see the `Quantum ESPRESSO official website <https://www.quantum-espresso.org/>`_.

Installing ELSES
---------------------

ELSES is a large-scale electronic structure calculation program, used in the fullerene MD tutorial (ELSES version).

The ELSES source code is provided to members of the `ELSES Consortium <http://www.elses.jp/>`_.
Participation in the consortium is required to use it.

For more details, see the `ELSES official website <http://www.elses.jp/>`_.

Installing odatse-aenet
----------------------------

.. code-block:: bash

   git clone <repository-url>
   cd odatse-aenet
   pip3 install .

After installation, the ``odatse-aenet`` command becomes available.

Verification
------------

.. code-block:: bash

   odatse-aenet --version

If the version number is displayed, the installation was successful.
