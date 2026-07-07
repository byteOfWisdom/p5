#!/usr/bin/python3
import std
import propeller as p
import numpy as np


def load_data(file):
    data = std.load_csv(file, "\t", 1)
    return (data[0], p.ev(data[1], np.sqrt(data[1])))


def line_func(x, mu, amp, sigma, a):
    return std.gaussian(x, amp, mu, sigma) + a


def get_peaks(x, y, guesses, plot=False):
    if plot:
        std.default.plt_errorbar(x, y)

    params, err, rsq = [], [], []
    for peak in guesses:
        p0 = [peak, ~y[peak], 50, 0]
        std.default.plt_errorbar(x, y)
        std.default.plt_errorbar(x[peak - 100: peak + 100], y[peak - 100: peak + 100])
        params_single, (err_single, rsq_single) = std.curve_fit(line_func, x[peak - 100: peak + 100], y[peak - 100: peak + 100], p0)
        std.default.plt_finish("", "")
        if plot:
          std.default.plt_func(line_func, params_single, xrange=(peak - 200, peak + 200))
        params.append(params_single)
        err.append(err_single)
        rsq.append(rsq_single)

    if plot:
        std.default.plt_finish("Bin", "Anzahl")
    return p.ev(params, err), rsq


if __name__ == "__main__":
    bin, count = load_data("../data/links/ba_ganz_lang.txt")
    get_peaks(bin, count, [420, 750, 930, 1190, 1670, 3870, 4230, 4970, 5390], True)
