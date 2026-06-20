# %%
import std
import tomllib
from matplotlib import pyplot as plt
import numpy as np
import propeller as p


# %%
fhandle = open("../data/meta.toml", "rb")
meta = tomllib.load(fhandle)
fhandle.close()


def quadratic(x, a, b, c):
    denom = (x ** 2) * b
    return (a / denom) + c

# %%
fname = "../data/" + meta["abstand_sq"]["messung_1"]["csv_file"]
data = std.util.load_csv(fname, skiprows=1)

distance = p.ev(data[3] + 11, 0.5)
dose_before = p.ev(data[1], 0.1)
dose_after = p.ev(data[2], 0.1)
time = p.ev(data[0], 0.5)

dose_rate = (dose_after - dose_before) / time

fit, (errs, rsq) = std.curve_fit(quadratic, distance, dose_rate)

plt.errorbar(~distance, ~dose_rate, xerr=p.ve(distance)[1], yerr=p.ve(dose_rate)[1], **std.default.error_bar_def)
xrange = np.linspace(min(~distance) - 3, max(~distance) + 3)
plt.plot(xrange, quadratic(xrange, *fit))
std.default.plt_finish("distance / cm", "dose rate / mSv/s")
