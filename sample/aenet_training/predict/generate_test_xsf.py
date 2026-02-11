"""Generate XSF test files for AENET prediction.

Creates N2 dimer structures at evenly spaced N-N distances
for evaluating the trained ANN potential.
"""

import os


def generate_test_xsf(output_dir="../predict_data_set_test",
                       d_min=0.0, d_max=2.0, d_step=0.01):
    """Generate XSF files with N2 dimers at various bond lengths.

    Parameters
    ----------
    output_dir : str
        Directory for output XSF files.
    d_min : float
        Minimum N-N distance (Angstrom).
    d_max : float
        Maximum N-N distance (Angstrom).
    d_step : float
        Distance step (Angstrom).
    """
    os.makedirs(output_dir, exist_ok=True)

    d = d_min
    while d <= d_max + 1e-9:
        filename = os.path.join(output_dir, f"test_{d:.2f}.xsf")
        with open(filename, "w") as f:
            f.write("ATOMS\n")
            f.write("N             0.0000000000        0.0000000000        0.0000000000\n")
            f.write(f"N             0.0000000000        0.0000000000        {d}\n")
        d += d_step


if __name__ == "__main__":
    generate_test_xsf()
