# %%
import propeller as p
import numpy as np
import std
from matplotlib import pyplot as plt
import scipy

# %%
def area_fraction(dist, detector_size):
    denominator = 2 * np.sqrt(1 + (detector_size / dist) ** 2)
    # denominator = 2 * np.sqrt(1 + (dist / detector_size) ** 2)
    return 0.5 - (1 / denominator)


def activity_at_time(starting_activity, time, half_life):
    mean_life =  half_life / np.log(2)
    return starting_activity * np.exp(- time / mean_life)


def time_with_err(time, err, unit):
    conversion = np.timedelta64(1, "D").item().total_seconds() * 365 if unit=="Y" else 1
    value =  time * conversion
    err = err * conversion
    return p.ev(value, err)


# %%
start_date = np.datetime64("2021-10-01")
experiment_date = np.datetime64("2026-04-29")
time_difference = experiment_date - start_date
delta_t = p.ev(time_difference.item().total_seconds(), np.timedelta64(1, "D").item().total_seconds())

# source: https://atom.kaeri.re.kr/nuchart/#
eu_152_hlt = time_with_err(13.517, 0.006, "Y")
cs_137_hlt = time_with_err(30.04, 0.04, "Y")
co_60_hlt =  time_with_err(5.2714, 0.0006, "Y")

activity = {
    "eu": activity_at_time(709e3, delta_t, eu_152_hlt),
    "cs": activity_at_time(405e3, delta_t, cs_137_hlt),
    "co": activity_at_time(67e3, delta_t, co_60_hlt),
}

print("152 Eu:", activity["eu"].format(), "Bq")
print("137 Cs:", activity["cs"].format(), "Bq")
print("60 Co:", activity["co"].format(), "Bq")

# %%
radius = {
    "scint": 50.8e-3 / 2,
    "ge": 55.7e-3 / 2,
}

duration = {
    "scint": 400,
    "ge": 300
}

dist = {
    "ge": {
        "eu": p.ev(14e-2, 5e-3),
        "cs": p.ev(85e-3, 5e-3),
        "co": p.ev(2e-3, 5e-3),
    },
    "scint": {
        "eu":  p.ev(170e-3, 5e-3),
        "cs":  p.ev(90e-3, 5e-3),
        "co":  p.ev(10e-3, 5e-3),
    }
}

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


for d, e in std.mesh(["ge", "scint"], ["cs", "co"]):
    counts_in_peak = fitted_peaks[d][e][0][0]
    transition_chance = p.ev(0.947, 0.002) # just looking at the ceasium peak here
    covered_area = area_fraction(dist[d][e], radius[d])
    emitted_rays = activity[e] * duration[d]
    reaching_detector = covered_area * emitted_rays * transition_chance
    ratio = counts_in_peak / reaching_detector

    print("running for:", d, e)
    print("area fraction is", covered_area.format())
    print("counts in peak are:", counts_in_peak.format())
    print("activity is:", activity[e].format())
    print("efficiency is:", ~ratio * 100, "%")
    print()

# %%
europium_line_matches = std.load_csv("../figs/europium_lit.csv", skiprows=1)
energy = europium_line_matches[0]
covered_area = area_fraction(dist["ge"]["eu"], radius["ge"])
total_gammas = activity["eu"] * duration["ge"]
intensity = europium_line_matches[1] / 100
reaching_detector = total_gammas * covered_area * intensity
ids = list(map(int, europium_line_matches[2]))
peak_areas = fitted_peaks["ge"]["eu"][0][ids]
effs = peak_areas / reaching_detector
# res, cov = scipy.optimize.curve_fit(lambda x, a, b, c, d: a * np.exp(- b * x) + c * np.exp(- d * x), ~energy, ~effs, p0=[1,1,1,1])
res, (err, rsq) = std.odr_fit(lambda x, a, b, c, d: a * np.exp(- b * x) + c * np.exp(- d * x), energy, effs)

print(effs)
plt.errorbar(~energy, ~effs, xerr=p.error(energy), yerr=p.error(effs), **std.default.error_bar_def)
std.default.plt_finish("Energie / keV", "$\\epsilon$")
