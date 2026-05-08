# %%
import propeller as p
import numpy as np
import std
from matplotlib import pyplot as plt


# %%
def area_fraction(dist, detector_size):
    denominator = 2 * np.sqrt(1 + (detector_size / dist) ** 2)
    return 1 / denominator


def activity_at_time(starting_activity, time, half_life):
    decay_const = np.log(2) / half_life
    return starting_activity * np.exp(- time * decay_const)


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
        "eu":  p.ev(10e-3, 5e-3),
        "cs":  p.ev(90e-3, 5e-3),
        "co":  p.ev(170e-3, 5e-3),
    }
}

elements = activity.keys()
detectors = radius.keys()

for d, e in std.mesh(detectors, elements):
    counts_in_peak = 1 # TODO: read correct data
    transition_chance = 1
    covered_area = area_fraction(dist[d][e], radius[d])
    emitted_rays = activity[e] * duration[d]
    reaching_detector = covered_area * emitted_rays * transition_chance
    ratio = counts_in_peak / reaching_detector
    print(ratio.format())
