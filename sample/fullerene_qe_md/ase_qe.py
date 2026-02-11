"""Fullerene NVT molecular dynamics using Quantum ESPRESSO via ASE."""

from ase import units
from ase.build import make_supercell
from ase.calculators.espresso import Espresso
from ase.io import read
from ase.io import Trajectory
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

# 1. Read fullerene structure from CIF
atoms = read("C.cif")
atoms.pbc = [1, 1, 1]
supercell_matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
supercell = make_supercell(atoms, supercell_matrix)

# 2. Quantum ESPRESSO settings
input_data = {
    'control': {
        'calculation': 'scf',
        'prefix': 'diamond_c',
        'pseudo_dir': './',
        'outdir': './',
    },
    'system': {
        'ecutwfc': 40,
        'ecutrho': 160,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.02,
    },
    'electrons': {
        'conv_thr': 1e-7,
    },
}

pseudopotentials = {'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF'}
command = 'srun -n 64 pw.x <espresso.pwi> espresso.pwo'

calc = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    kpts=(1, 1, 1),
    tprnfor=True,
    input_data=input_data,
)
supercell.calc = calc

# 3. MD settings
dt = 2 * units.fs
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
