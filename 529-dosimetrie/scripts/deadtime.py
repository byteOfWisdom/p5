import std
import propeller as p
import tomllib
import numpy as np

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["totzeit"]
fhandle.close()


def non_paralysing(x, a, b, tau):
    m = tau + (1 / (a * x + b))
    return 1 / m
    # return a * x / (1 - a * tau * x)


def paralysing(x, a, b):
    return x


print(metadata)

data = std.util.load_csv("../data/" + metadata["csv_file"], skiprows=1)
emission_current = data[0]
countrate = p.ev(data[1], np.sqrt(data[1]))
std.default.plt_errorbar(emission_current, countrate, "gemessene Zählraten")
std.default.plt_finish("Heizstrom $I_C$ / mA", "Ionendosisleistung / $Akg^{-1}$")


countrate = countrate#[emission_current <= 0.2]
emission_current = emission_current#[emission_current <= 0.2]
# dose_rate = current_to_countrate(emission_current)#[emission_current <= 0.2]

res, (err, rsq) = std.curve_fit(non_paralysing, emission_current, countrate)


std.default.plt_errorbar(emission_current, countrate)
std.default.plt_func(non_paralysing, res, xrange=(0.01, None))
std.default.plt_finish("Heizstrom $I_C$ / mA", "Ionendosisleistung / $Akg^{-1}$")
