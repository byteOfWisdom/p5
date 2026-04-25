# from scipy import curve_fit
import numpy as np
import std
from matplotlib import pyplot as plt

def fit_biggest_peak(channel, amplitude):
    max_loc = channel[np.where(amplitude == np.max(amplitude))]
    p0 = [amplitude[max_loc], max_loc, 1]
    # fit_res = curve_fit(channel, amplitude, std.gaussian, p0=p0)
    fit_res, _ = std.fit_func(std.gaussian, channel, amplitude, p0=p0, force_cf=True)
    return fit_res


def decomp_spectrum(channel, amplitude):
    plt.plot(channel, amplitude)
    res = fit_biggest_peak(channel, amplitude)
    plt.plot(channel, std.gaussian(channel, *res), linestyle="dashed")
    plt.show()
