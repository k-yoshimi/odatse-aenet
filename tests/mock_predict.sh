#!/bin/sh
# Mock predict.x for testing
# Reads predict.in and structure_data.xsf, outputs a fake predict.out
cat > predict.out << 'PREDICT_OUTPUT'

 ======================================================================
                predict.x - energy/force prediction
 ======================================================================

 Number of atoms   :         2
 Number of species :         1

 Species : N

 Atomic Energy Network Potentials
 --------------------------------

 N  (  5- 5) N.5t-5t.ann

 Energy (Hartree)  :    -0.12345678
 Total energy      :    -3.36200000 eV
 Energy/atom       :    -1.68100000 eV

PREDICT_OUTPUT
