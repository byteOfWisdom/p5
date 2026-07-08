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
        "na": [511e3],
        'ba': [32.19e3, np.nan, np.nan, 81e3, np.nan, np.nan, np.nan, 356e3, 383.8e3]
    },
    "rechts": {
        # "na": [511e3, 1274.5e3, np.nan],
        "na": [511e3],
        'ba': [32.19e3, np.nan, np.nan, 81e3, np.nan, np.nan, np.nan, 356e3, 383.8e3]
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
    params, (err, rsq) = std.curve_fit(std.linear, energies, lines)
    std.default.plt_errorbar(energies, lines)
    std.default.plt_func(std.linear, params, f"R^2 = {round(rsq, 4)}")
    std.default.plt_finish("Energie / keV", "$\\mu$ / Bins")
    table = {
        "$\\mu$ / Bins": lines,
        "Energie / keV": energies,
    }
    std.print_tex_table(table, f"../latex/{side}_energy_cal.table")


if __name__ == "__main__":
    run_side("links")
    run_side("rechts")
