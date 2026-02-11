"""Plot kinetic energy per atom from ASE MD log file."""

import argparse

import matplotlib.pyplot as plt

BOLTZMANN_EV = 8.617333e-5  # Boltzmann constant [eV/K]


def parse_md_log(filepath, natoms):
    """Parse an ASE MD log file and extract kinetic energy per atom.

    Parameters
    ----------
    filepath : str
        Path to the ASE MD log file (e.g., ``md.log``).
    natoms : int
        Number of atoms in the system, used to normalize kinetic energy.

    Returns
    -------
    time : list of float
        Simulation time in femtoseconds.
    ekin : list of float
        Kinetic energy per atom in eV/atom.
    """
    time = []
    ekin = []
    with open(filepath, "r") as f:
        for line in f:
            if "Time[ps]" in line or line.strip() == "":
                continue
            cols = line.split()
            time.append(float(cols[0]) * 1000)  # ps -> fs
            ekin.append(float(cols[3]) / natoms)
    return time, ekin


def plot_ekin(time, ekin, temperature, output_file):
    """Plot kinetic energy per atom with a Boltzmann equipartition reference.

    Draws the kinetic energy time series and a horizontal reference line
    at ``(3/2) * k_B * T``.

    Parameters
    ----------
    time : list of float
        Simulation time in femtoseconds.
    ekin : list of float
        Kinetic energy per atom in eV/atom.
    temperature : float
        Target temperature in Kelvin for the reference line.
    output_file : str
        Path to save the plot image.
    """
    ekin_ref = (3 / 2) * BOLTZMANN_EV * temperature

    plt.figure(figsize=(8, 6))
    plt.plot(time, ekin, marker='o', linestyle='-', color='b', label='Ekin')
    plt.axhline(ekin_ref, color='red', linestyle='--',
                label=f"Ekin (T={temperature} K): {ekin_ref:.4f} eV/atom")
    plt.xlabel("Time [fs]", fontsize=14)
    plt.ylabel("Ekin [eV/atom]", fontsize=14)
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Ekin from MD log")
    parser.add_argument("--log", default="md.log", help="MD log file path")
    parser.add_argument("--natoms", type=int, default=480, help="Number of atoms")
    parser.add_argument("--temperature", type=float, default=300, help="Temperature [K]")
    parser.add_argument("--output", default="time_energy.png", help="Output image file")
    args = parser.parse_args()

    time, ekin = parse_md_log(args.log, args.natoms)
    plot_ekin(time, ekin, args.temperature, args.output)
