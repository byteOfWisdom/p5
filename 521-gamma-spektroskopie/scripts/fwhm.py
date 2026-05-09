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
for detector, element in std.mesh(["ge", "scint"], ["eu", "cs", "co"]):
    fwhm = fitted_peaks[detector][element][sigma] * fwhm_const
    fwhm_energy = fwhm * energy_cal[detector]["slope"]
    print(f"{detector}, {element}, {fwhm}, {fwhm_energy}")


# %%
def func(e, c, a):
    return a + c * (e ** 0.5)

# %%
fwhm = fitted_peaks["ge"]["eu"][sigma] * fwhm_const
fwhm_energy = fwhm * energy_cal["ge"]["slope"]
energy = fitted_peaks["ge"]["eu"][mu] * energy_cal["ge"]["slope"] + energy_cal["ge"]["offset"]

europium_line_matches = std.load_csv("../figs/europium_lit.csv", skiprows=1)
ids = list(map(int, europium_line_matches[2]))
fwhm_energy = fwhm_energy[ids]
energy = energy[ids]
# fwhm = fitted_peaks["scint"]["eu"][sigma] * fwhm_const
# fwhm_energy = fwhm * energy_cal["scint"]["slope"]
# energy = fitted_peaks["scint"]["eu"][mu] * energy_cal["scint"]["slope"] + energy_cal["scint"]["offset"]

res, (err, rsq) = std.fit_func(func, ~energy, ~fwhm_energy, force_cf=True)
# res, (err, rsq) = std.fit_func(lambda x, a, b:a * x + b, ~energy, ~fwhm_energy, force_cf=True)
print(rsq)
# res, _ = scipy.optimize.curve_fit(func, ~energy, ~fwhm_energy)

plt.errorbar(~energy, ~fwhm_energy, p.error(fwhm_energy), p.error(energy), **std.default.error_bar_def)
# plt.scatter(~energy, ~fwhm_energy)
erange = np.linspace(min(~energy), max(~energy), 10000)
plt.plot(erange, func(erange, *res), label="fit")
std.default.plt_finish("Energie / eV", "FWHM / eV")
