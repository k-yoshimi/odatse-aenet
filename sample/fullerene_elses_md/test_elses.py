"""Fullerene NVT molecular dynamics using ELSES via ASE."""

from ase import units
from ase.build import make_supercell
from ase.io import read
from ase.io import Trajectory
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from ElsesCalculator import elses

# 1. Read fullerene structure from CIF
atoms = read("C.cif")
atoms.pbc = [1, 1, 1]
supercell_matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
supercell = make_supercell(atoms, supercell_matrix)

print("Number of atoms:", len(supercell))
print("PBC:", supercell.get_pbc())

# 2. ELSES calculator setup (adjust command for your environment)
command = 'elses-xml-generate generate.xml C480.xml ; srun elses config.xml > log.txt'
calc = elses(command=command)
supercell.calc = calc

# 3. MD settings
dt = 1.0 * units.fs
temperature_K = 300
nsteps = 100

MaxwellBoltzmannDistribution(supercell, temperature_K=temperature_K)

taut = 0.5 * 10 * units.fs
traj = Trajectory("benzene_optimization.traj", "w", supercell)
dyn = NVTBerendsen(
    supercell, timestep=dt, temperature=temperature_K,
    taut=taut, trajectory=traj, logfile='md.log',
)

dyn.run(nsteps)
