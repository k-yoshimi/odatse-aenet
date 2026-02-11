import xml.etree.ElementTree as ET

import numpy as np
from ase.calculators.calculator import FileIOCalculator, PropertyNotPresent
from ase.units import Hartree, Bohr
from collections import Counter

# Unit conversion constants (from ase.units)
HARTREE_TO_EV = Hartree              # 27.2114 eV
BOHR_TO_ANGSTROM = Bohr              # 0.52918 Angstrom
AU_FORCE_TO_EV_ANGSTROM = Hartree / Bohr

# Keep old names for backward compatibility
ev = HARTREE_TO_EV
angst = BOHR_TO_ANGSTROM
eV_angst = AU_FORCE_TO_EV_ANGSTROM


class ElsesCalculator(FileIOCalculator):
    """ASE FileIOCalculator for the ELSES electronic structure code.

    Generates XYZ input for ELSES, executes the external command, and
    reads back energy from ``Output.txt`` and forces from ``restart.xml``.

    Parameters
    ----------
    restart : str, optional
        Path to restart file.
    label : str, optional
        Prefix for input/output files (default: ``'elses'``).
    atoms : ase.Atoms, optional
        ASE Atoms object.
    **kwargs
        Additional keyword arguments passed to
        :class:`ase.calculators.calculator.FileIOCalculator`.
    """

    implemented_properties = ['energy', 'forces']

    def __init__(self, restart=None, label='elses', atoms=None, **kwargs):
        super().__init__(restart=restart, label=label, atoms=atoms, **kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        """Generate XYZ input file for ELSES from an ASE Atoms object.

        Creates a file named ``<formula>.xyz`` containing the number of atoms,
        cell parameters, and atomic positions.

        Parameters
        ----------
        atoms : ase.Atoms
            Atoms object with positions and cell information.
        properties : list of str, optional
            Requested properties (unused, kept for API compatibility).
        system_changes : list of str, optional
            System changes since last call (unused, kept for API compatibility).
        """
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        cell_params = atoms.cell.cellpar()

        element_counts = Counter(symbols)
        formula = ''.join(
            f"{el}{n}" if n > 1 else el
            for el, n in element_counts.items()
        )

        with open(f"{formula}.xyz", "w") as f:
            f.write(f"{len(symbols)}  "
                    f"{cell_params[0]:.1f} {cell_params[1]:.1f} {cell_params[2]:.1f}\n")
            f.write(f"{formula}\n")
            for symbol, pos in zip(symbols, positions):
                f.write(f"{symbol:<4} {pos[0]: .15E} {pos[1]: .15E} {pos[2]: .15E}\n")

    def read_results(self):
        """Read energy from ``Output.txt`` and forces from ``restart.xml``.

        Forces are converted from atomic units to eV/Angstrom.
        Energy per atom is read from the last ``Energy summary`` line in
        ``Output.txt`` and multiplied by the number of atoms.

        Returns
        -------
        energy : float
            Total energy in eV.
        forces : numpy.ndarray
            Forces array of shape ``(natoms, 3)`` in eV/Angstrom.

        Raises
        ------
        ase.calculators.calculator.PropertyNotPresent
            If no energy summary is found in ``Output.txt``.
        """
        # Force extraction from restart.xml
        tree = ET.parse('restart.xml')
        root = tree.getroot()
        total_atoms = int(root.find('.//total_number_of_atoms').text)

        forces_au = []
        for atom in root.findall('atom'):
            force_text = atom.find('force').text.strip()
            forces_au.append([float(x) for x in force_text.split()])
        self.results['forces'] = np.array(forces_au) * AU_FORCE_TO_EV_ANGSTROM

        # Energy extraction from Output.txt
        energy_values = []
        with open('Output.txt', 'r') as f:
            for line in f:
                if 'Energy summary (eV/atom):' in line:
                    energy_values.append(float(line.split()[4]))

        if not energy_values:
            raise PropertyNotPresent(
                'Property "Energy" not available. Please try running\n'
                'ELSES first by calling Atoms.get_potential_energy().'
            )

        self.results['energy'] = energy_values[-1] * total_atoms
        return self.results['energy'], self.results['forces']


# Backward-compatible alias
elses = ElsesCalculator
