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
        "ug": std.load_csv("../data/undergrd_ge.txt")
    },
    "scint": {
        "eu": std.load_csv("../data/eu_nai.txt"),
        "cs": std.load_csv("../data/cs_nai_fixed.txt"),
        "co": std.load_csv("../data/co_nai.txt"),
        "ug": std.load_csv("../data/undergrd_nai.txt")
    }
}

# %%
for detector, element in std.mesh(["ge", "scint"], ["cs", "co"]):
    print(detector, element)
    total_hits = sum(np.vectorize(int)(spectra[detector][element][1]) - np.vectorize(int)(spectra[detector]["ug"][1]))
    print(total_hits)
    peak_to_total = fitted_peaks[detector][element][0] / total_hits
    print("mu:", fitted_peaks[detector][element][1])
    print("PTT:", peak_to_total)
    print(sum(peak_to_total).format())
    print()


# %%
print()
