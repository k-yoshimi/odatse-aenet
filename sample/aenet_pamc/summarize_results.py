"""Summarize PAMC weight results for each temperature step.

Reads ``weight.txt`` files from ODAT-SE PAMC output directories and
produces per-temperature summary files with both normalized and
non-normalized weights.
"""

import argparse
import datetime
import os
import time

import numpy as np

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


def load_weights(output_dirname, num_procs, rep_per_pro, Tnum, n_search_val, idnum):
    """Load weight data from PAMC output and organize by temperature step.

    Parameters
    ----------
    output_dirname : str
        Path to the ODAT-SE output directory.
    num_procs : int
        Number of MPI processes used in the PAMC run.
    rep_per_pro : int
        Number of replicas per process.
    Tnum : int
        Number of temperature steps.
    n_search_val : int
        Number of search variables (dimension).
    idnum : int
        Target data ID number to extract.

    Returns
    -------
    T : numpy.ndarray
        Temperature values for each step, shape ``(Tnum,)``.
    txt_datas : numpy.ndarray
        Aggregated data array, shape ``(Tnum, nreplica, outfile_column)``.
    """
    outfile_column = 5 + n_search_val
    infile_column = 6 + n_search_val
    nreplica = num_procs * rep_per_pro

    txt_datas = np.zeros([Tnum, nreplica, outfile_column])
    counter_eachT = np.zeros(Tnum, dtype=int)
    T = -np.ones(Tnum, dtype=np.float64)

    ind_bool = np.ones(infile_column, dtype=bool)
    ind_bool[1] = False

    for i in range(num_procs):
        f = np.loadtxt(f"{output_dirname}/{i}/weight.txt")
        f_in_specific_idnum = f[:, 3] == idnum
        temp_txt_datas_eachT = f[f_in_specific_idnum, :][:, ind_bool]
        temp_txt_datas_T = f[f_in_specific_idnum, :][:, 1]

        for j in range(Tnum):
            if counter_eachT[j] == 0:
                T[j] = temp_txt_datas_T[temp_txt_datas_eachT[:, 0] == j][0]
            temp_txt_datas = temp_txt_datas_eachT[temp_txt_datas_eachT[:, 0] == j, 1:]
            temp_nrep = temp_txt_datas.shape[0]
            txt_datas[j, counter_eachT[j]:counter_eachT[j] + temp_nrep, 0] = i
            txt_datas[j, counter_eachT[j]:counter_eachT[j] + temp_nrep, 1:] = temp_txt_datas
            counter_eachT[j] += temp_nrep

    return T, txt_datas


def save_summary(txt_datas, T, Tnum, n_search_val, out_dir, dataname, idnum, date, fmt,
                 normalize=False):
    """Save per-temperature summary files.

    Parameters
    ----------
    txt_datas : numpy.ndarray
        Aggregated data array from :func:`load_weights`.
    T : numpy.ndarray
        Temperature array.
    Tnum : int
        Number of temperature steps.
    n_search_val : int
        Number of search variables.
    out_dir : str
        Output directory path.
    dataname : str
        Dataset name for filenames.
    idnum : int
        Target data ID number.
    date : str
        Date string for headers.
    fmt : str
        Format string for ``numpy.savetxt``.
    normalize : bool, optional
        Whether to normalize weights (default: False).
    """
    os.makedirs(out_dir, exist_ok=True)

    for i in range(Tnum):
        if normalize:
            sum_wgt = np.sum(txt_datas[i, :, 4])
            txt_datas[i, :, 4] = txt_datas[i, :, 4] / sum_wgt
            norm_note = (
                f"Weights are normalized. ( sum(weight) = {np.sum(txt_datas[i, :, 4])} )\n"
                f"sum(non_normalized_weight) = {sum_wgt} \n"
            )
        else:
            norm_note = "Weights are not normalized.\n"

        header = (
            f"generate date : {date}\n"
            f"T = {T[i]} \n"
            f"{norm_note}\n"
            f"#0 n_proc\n#1 n_replica\n#2 id_num\n#3 fx\n#4 weight\n"
        )
        for j in range(n_search_val):
            header += f"#{j + 5} value_{j + 1}\n"

        suffix = "norm_wgt" if normalize else "no_norm_wgt"
        filename = f"{i:03}_Tstep{i}_{dataname}_idnum{idnum}_Tnum{Tnum}_fx_{suffix}.txt"
        np.savetxt(os.path.join(out_dir, filename), txt_datas[i], header=header, fmt=fmt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarize PAMC weight results for each temperature step.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-d", "--dataname", help="Dataset name for output filenames")
    parser.add_argument("-i", "--input_filename", help="input.toml used for ODAT-SE", required=True)
    parser.add_argument(
        "-t", "--specific_Tnum",
        help="Process only this many temperature steps (default: from input.toml)",
        default=None,
    )
    parser.add_argument("--idnum", help="Target data ID number", required=True, type=int)
    args = parser.parse_args()

    date = datetime.date.today().strftime("%Y%m%d")
    with open(args.input_filename, "rb") as f:
        toml_dir = tomllib.load(f)

    output_dirname = toml_dir["base"].get("output_dir", ".")
    n_search_val = toml_dir["base"]["dimension"]
    rep_per_pro = toml_dir["algorithm"]["pamc"]["replica_per_prop"]
    Tnum = int(args.specific_Tnum) if args.specific_Tnum else toml_dir["algorithm"]["pamc"]["Tnum"]
    num_procs = len([
        f for f in os.listdir(output_dirname)
        if os.path.isdir(os.path.join(output_dirname, f))
    ])

    fmt = "%i %i %i %.15e %.15e" + " %.15e" * n_search_val

    T, txt_datas = load_weights(
        output_dirname, num_procs, rep_per_pro, Tnum, n_search_val, args.idnum,
    )

    base = f"{date}_{args.dataname}_idnum{args.idnum}_Tnum{Tnum}"
    save_summary(
        txt_datas.copy(), T, Tnum, n_search_val,
        f"{base}_fx_non_norm_wgt_eachT", args.dataname, args.idnum, date, fmt,
        normalize=False,
    )
    save_summary(
        txt_datas, T, Tnum, n_search_val,
        f"{base}_fx_normwgt_eachT", args.dataname, args.idnum, date, fmt,
        normalize=True,
    )
