"""Generate QE input files for N2 dimer with random bond lengths."""

import os
import random
import shutil


def generate_structures(
    n_structures=20,
    r_min=0.5,
    r_max=2.0,
    template_file="template.txt",
    pseudo_dir="./pseudo_Potential/",
    output_prefix="directory",
):
    """Generate Quantum ESPRESSO input directories with random N-N distances.

    For each structure, a directory is created by copying the pseudopotential
    directory, and the template file is filled with a random N-N bond length.

    Parameters
    ----------
    n_structures : int, optional
        Number of structures to generate (default: 20).
    r_min : float, optional
        Minimum N-N bond length in Angstrom (default: 0.5).
    r_max : float, optional
        Maximum N-N bond length in Angstrom (default: 2.0).
    template_file : str, optional
        Path to the template input file containing ``'value_01'`` placeholder
        (default: ``'template.txt'``).
    pseudo_dir : str, optional
        Directory containing pseudopotential files to copy
        (default: ``'./pseudo_Potential/'``).
    output_prefix : str, optional
        Prefix for output directory names (default: ``'directory'``).
    """
    placeholders = ["value_01"]
    distances = [random.uniform(r_min, r_max) for _ in range(n_structures)]

    for i, distance in enumerate(distances):
        dir_name = f"{output_prefix}_{i}"
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
        shutil.copytree(pseudo_dir, dir_name)

        with open(template_file, "r") as f:
            content = f.read()
        for placeholder, value in zip(placeholders, [distance]):
            content = content.replace(placeholder, str(value))

        output_path = os.path.join(dir_name, "n2_dimer.pwi")
        with open(output_path, "w") as f:
            f.write(content)


if __name__ == "__main__":
    generate_structures()
