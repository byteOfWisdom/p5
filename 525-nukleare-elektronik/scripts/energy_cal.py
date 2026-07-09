#!/usr/bin/python3
import std
import propeller as p
import numpy as np


files = {
    "rechts": ["na", "ba"],
    "links": ["na", "ba"],
}

energy = {
    "links": {
        # "na": [511e3, 1274.5e3],
        "na": [511],
        'ba': [32.19, np.nan, np.nan, 81, np.nan, np.nan, np.nan, 356, 383.8]
    },
    "rechts": {
        # "na": [511e3, 1274.5e3, np.nan],
        "na": [511],
        'ba': [32.19, np.nan, np.nan, 81, np.nan, np.nan, np.nan, 356, 383.8]
    }
}


def load_data(side, element):
    fname = f"../data/fits/{side}/{element}_ganz.csv"
    data = std.load_csv(fname, skiprows=1)

    filter = np.array(energy[side][element]) == np.array(energy[side][element])
    return (data[1])[filter], np.array(energy[side][element])[filter]


def run_side(side):
    na_lines, na_energies = load_data(side, "na")
    ba_lines, ba_energies = load_data(side, "ba")
    lines = np.append(na_lines, ba_lines)
    energies = np.append(na_energies, ba_energies)
    params, (err, rsq) = std.odr_fit(std.linear, lines, energies)
    res = p.ev(params, err)
    std.write_file(f"../data/fits/{side}/energy_cal_res.txt", f"{std.si_string("a", res[0], "\\kilo\\eV\\per\\bin")} und {std.si_string("b", res[1], "\\kilo\\eV")}")
    std.default.plt_errorbar(lines, energies, marker="x")
    std.default.plt_func(std.linear, params, f"R^2 = {round(rsq, 4)}")
    std.default.plt_finish("$\\mu$ / Bins", "Energie / keV")
    table = {
        "$\\mu$ / Bins": lines,
        "Energie / keV": energies,
    }
    std.print_tex_table(table, f"../latex/{side}_energy_cal.table")
    print(p.ev(params, err))


if __name__ == "__main__":
    run_side("links")
    run_side("rechts")
