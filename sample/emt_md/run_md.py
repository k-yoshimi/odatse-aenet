"""Simple NVT molecular dynamics of Ni bulk using ASE's EMT calculator.

Demonstrates a two-stage temperature ramp: equilibration at low temperature
followed by a production run at a higher temperature.
"""

import argparse

import numpy as np
from ase import units
from ase.build import bulk, make_supercell
from ase.calculators.emt import EMT
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


def run_md(element="Ni", supercell_size=3, temp_equil=10, nsteps_equil=200,
           temp_prod=500, nsteps_prod=400, dt_fs=2.0, taut_fs=20.0,
           traj_file="md.traj"):
    """Run two-stage NVT MD simulation with EMT potential.

    Parameters
    ----------
    element : str, optional
        Element symbol (default: ``'Ni'``).
    supercell_size : int, optional
        Supercell multiplier in each direction (default: 3).
    temp_equil : float, optional
        Equilibration temperature in K (default: 10).
    nsteps_equil : int, optional
        Number of equilibration steps (default: 200).
    temp_prod : float, optional
        Production temperature in K (default: 500).
    nsteps_prod : int, optional
        Number of production steps (default: 400).
    dt_fs : float, optional
        Time step in femtoseconds (default: 2.0).
    taut_fs : float, optional
        Berendsen thermostat coupling time in fs (default: 20.0).
    traj_file : str, optional
        Output trajectory file (default: ``'md.traj'``).
    """
    atoms = bulk(element, "fcc", a=3.5, cubic=True)
    atoms = make_supercell(atoms, np.diag([supercell_size] * 3))
    atoms.calc = EMT()

    dt = dt_fs * units.fs
    taut = taut_fs * units.fs

    MaxwellBoltzmannDistribution(atoms, temperature_K=temp_equil)
    dyn = NVTBerendsen(atoms, dt, temp_equil, taut=taut, trajectory=traj_file)

    def print_status():
        print(f"time={dyn.get_time() / units.fs: 5.0f} fs "
              f"T={atoms.get_temperature(): 3.0f} K")

    dyn.attach(print_status, interval=20)

    print(f"--- Equilibration at {temp_equil} K ---")
    dyn.run(nsteps_equil)

    print(f"--- Production at {temp_prod} K ---")
    dyn.set_temperature(temp_prod)
    dyn.run(nsteps_prod)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ni bulk NVT MD with EMT calculator")
    parser.add_argument("--element", default="Ni", help="Element symbol")
    parser.add_argument("--supercell", type=int, default=3, help="Supercell multiplier")
    parser.add_argument("--temp-equil", type=float, default=10, help="Equilibration temperature [K]")
    parser.add_argument("--nsteps-equil", type=int, default=200, help="Equilibration steps")
    parser.add_argument("--temp-prod", type=float, default=500, help="Production temperature [K]")
    parser.add_argument("--nsteps-prod", type=int, default=400, help="Production steps")
    parser.add_argument("--traj", default="md.traj", help="Output trajectory file")
    args = parser.parse_args()

    run_md(
        element=args.element, supercell_size=args.supercell,
        temp_equil=args.temp_equil, nsteps_equil=args.nsteps_equil,
        temp_prod=args.temp_prod, nsteps_prod=args.nsteps_prod,
        traj_file=args.traj,
    )
