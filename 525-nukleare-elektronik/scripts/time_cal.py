import std
import numpy as np
import propeller as p


def load_data():
    fname = "../data/fits/prompt.csv"
    data = std.load_csv(fname, skiprows=1)
    return data[1]

def time_cal(plot=False):
    lines = load_data()
    times = 16 * np.arange(5) + 8
    res, (err, rsq) = std.curve_fit(std.linear, times, lines)
    # res, (err, rsq) = std.odr_fit(std.linear, lines,  times)
    if plot:
        std.default.plt_errorbar(times, lines)
        # std.default.plt_errorbar(lines, times)
        std.default.plt_func(std.linear, res, f"$R^2 = {round(rsq, 4)}$")
        std.default.plt_finish("Zeit / ns", "Linie / Bins")
        # std.default.plt_finish("Linie / Bins", "Zeit / ns")
    return p.ev(res, err)


def invert(a, b):
    return 1 / a, -b / a


if __name__ == "__main__":
    print(time_cal(True))
