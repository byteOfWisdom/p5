# %%
import std
import numpy as np
from matplotlib import pyplot as plt
import propeller as p
import scipy


# %%
energy_cal = {
    "ge": {
        "slope": p.from_string("92.6664(19)"),
        "offset": p.from_string("154(19)")
    },
    "scint": {
        "slope": p.from_string("120.3(1.3)"),
        "offset": p.from_string("0.0(9.7)e+03")
    }
}

area, mu, sigma = 0, 1, 2

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

fwhm_const = np.sqrt(8 * np.log(2))

 # %%
for detector, element in std.mesh(["ge", "scint"], ["cs", "co"]):
    fwhm = fitted_peaks[detector][element][sigma] * fwhm_const
    fwhm_energy = fwhm * energy_cal[detector]["slope"]
    energy = fitted_peaks[detector][element][mu] * energy_cal[detector]["slope"] + energy_cal[detector]["offset"]
    for i in range(len(fwhm)):
        print(f"{detector} & {element} & {energy[i].format()} & {fitted_peaks[detector][element][sigma][i].format()} & {fwhm[i].format()} & {fwhm_energy[i].format()} \\\\")


# %%
def func(x, a, b):
    return a * x + b

# %%
fwhm = fitted_peaks["ge"]["eu"][sigma] * fwhm_const
fwhm_energy = fwhm * energy_cal["ge"]["slope"]
energy = fitted_peaks["ge"]["eu"][mu] * energy_cal["ge"]["slope"] + energy_cal["ge"]["offset"]

# fwhm = fitted_peaks["scint"]["eu"][sigma] * fwhm_const
# fwhm_energy = fwhm * energy_cal["scint"]["slope"]
# energy = fitted_peaks["scint"]["eu"][mu] * energy_cal["scint"]["slope"] + energy_cal["scint"]["offset"]

europium_line_matches = std.load_csv("../figs/europium_lit.csv", skiprows=1)
ids = list(map(int, europium_line_matches[2]))
fwhm_energy = fwhm_energy[ids]
fwhm = fwhm[ids]
energy = energy[ids]
energy = europium_line_matches[0] * 1e3

x = energy
y = fwhm_energy ** 2

res, (err, rsq) = std.curve_fit(func, x, y)
# res, (err, rsq) = std.fit_func(lambda x, a, b:a * x + b, ~energy, ~fwhm_energy, force_cf=True)
print(p.ev(res, err))
print(np.sqrt(p.ev(res, err)))
# res, _ = scipy.optimize.curve_fit(func, ~energy, ~fwhm_energy)

plt.errorbar(~x, ~y, p.error(y), p.error(x), **std.default.error_bar_def)
erange = np.linspace(0.8 * min(~energy), 1.05 * max(~energy), 10000)
plt.plot(erange, func(erange, *res), label=f"$R^2={round(rsq, 3)}$")
std.default.plt_finish("$E_\\gamma$ / eV", "$\\text{FWHM}^2$ / $\\text{eV}^2$", "../figs/fwhm_fit.pdf")

# %%
table = {
    "$E_\\gamma$ / eV": energy,
    "FWHM / Kanäle": fwhm,
    "FWHM / eV": fwhm_energy
}

print(table)
std.print_tex_table(table, "../figs/fwhm_eu_ge.table")
