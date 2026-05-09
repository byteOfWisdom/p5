# %%
import propeller as p
import numpy as np
import std
from matplotlib import pyplot as plt

fitted_peaks = {
    "ge": {
        "eu": std.load_csv("../figs/eu_ge_bin.csv", skiprows=1),
        "cs": std.load_csv("../figs/cs_ge_bin.csv", skiprows=1),
        "co": std.load_csv("../figs/co_ge_bin.csv", skiprows=1),
    },
    "scint": {
        "eu": std.load_csv("../figs/eu_nai_bin.csv", skiprows=1),
        "cs": std.load_csv("../figs/cs_nai_bin.csv", skiprows=1),
        "co": std.load_csv("../figs/co_nai_bin.csv", skiprows=1),
    }
}

spectra = {
    "ge": {
        "eu": std.load_csv("../data/eu_ge.txt"),
        "cs": std.load_csv("../data/cs_ge.txt"),
        "co": std.load_csv("../data/co_ge.txt"),
    },
    "scint": {
        "eu": std.load_csv("../data/eu_nai.txt"),
        "cs": std.load_csv("../data/cs_nai_fixed.txt"),
        "co": std.load_csv("../data/co_nai.txt"),
    }
}

# %%
for detector, element in std.mesh(["ge", "scint"], ["eu", "cs", "co"]):
    total_hits = sum(spectra[detector][element][1])
    peak_to_toal = fitted_peaks[detector][element][0] / total_hits
    print(peak_to_toal)
    print(sum(peak_to_toal).format())
