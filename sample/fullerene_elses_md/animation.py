"""Create animated GIF of MD trajectory with kinetic energy plot."""

import argparse

import numpy as np
from ase import units
from ase.io.trajectory import Trajectory
from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def create_animation(traj_file, output_file="animation.gif", natoms=480,
                     dt_fs=10.0, nframes=10):
    """Create an animated GIF showing structure and kinetic energy over time.

    Parameters
    ----------
    traj_file : str
        Path to the ASE trajectory file.
    output_file : str, optional
        Output GIF filename (default: ``'animation.gif'``).
    natoms : int, optional
        Number of atoms for kinetic energy normalization (default: 480).
    dt_fs : float, optional
        Time step in femtoseconds (default: 10.0).
    nframes : int, optional
        Number of frames in the animation (default: 10).
    """
    traj = Trajectory(traj_file)
    nsteps = len(traj) - 1

    time_array = np.arange(nsteps + 1) * dt_fs * units.fs
    ekin = [atoms.get_kinetic_energy() / natoms for atoms in traj]

    fig, ax = plt.subplots(1, 2, figsize=(9, 3), tight_layout=True)

    ekin_ref = 0.0388  # Ekin at T=300K [eV/atom]

    def update(iframe):
        idx = int(nsteps * iframe / nframes)

        ax[0].cla()
        ax[0].set_xlabel("time (fs)")
        ax[0].set_ylabel("Ekin (eV/atom)")
        ax[0].set_xlim(0, time_array[-1])
        ax[0].axhline(y=ekin_ref, color="red", linestyle="--", linewidth=1.5,
                       label=f"Ekin(T=300K): {ekin_ref} eV/atom")
        ax[0].legend(loc="lower right")
        ax[0].plot(time_array, ekin)
        ax[0].plot(time_array[idx], ekin[idx], marker="X", markersize=10)

        ax[1].cla()
        ax[1].set_title("Structure")
        ax[1].axis("off")
        plot_atoms(traj[idx], ax=ax[1], rotation="90x,0y", scale=0.5, radii=0.5)

        return ax[0].lines + ax[1].patches

    ani = FuncAnimation(fig, update, np.arange(nframes), blit=True, interval=250.0)
    ani.save(output_file, writer="imagemagick")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create MD trajectory animation")
    parser.add_argument("--traj", default="benzene_optimization.traj",
                        help="ASE trajectory file")
    parser.add_argument("--output", default="animation.gif", help="Output GIF file")
    parser.add_argument("--natoms", type=int, default=480, help="Number of atoms")
    parser.add_argument("--dt", type=float, default=10.0, help="Time step [fs]")
    parser.add_argument("--nframes", type=int, default=10, help="Number of animation frames")
    args = parser.parse_args()

    create_animation(args.traj, args.output, args.natoms, args.dt, args.nframes)
