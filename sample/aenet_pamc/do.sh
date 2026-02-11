#!/bin/sh

export PYTHONUNBUFFERED=1
export OMPI_MCA_rmaps_base_oversubscribe=true

#CMD="python3 ../../src/main.py"
CMD="odatse-aenet"

#MPIEXEC=""
MPIEXEC="mpiexec -np 2"

time $MPIEXEC $CMD input.toml
