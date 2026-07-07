import std
import numpy as np
import propeller as p
import time_cal
import scipy.special as sp
# import iminuit as im


def lt_spectrum(t, amplitude, tau, t0,  sigma):
    dt = t - t0
    # b = ((dt / sigma) - (sigma / tau)) / np.sqrt(2)
    b = (tau * dt - sigma**2) / (np.sqrt(2) * sigma * tau)
    # return b
    # d = (sigma ** 2 + tau * t0) / (np.sqrt(2) * sigma * tau)
    c = amplitude * (1 + sp.erf(b)) / (2 * tau)
    return c * np.exp(((sigma**2) - (2 * tau * dt)) / (2 * tau ** 2))


if __name__ == "__main__":
    time_cal_params = time_cal.invert(*time_cal.time_cal())
    time_slope = time_cal_params[0]
    data = std.load_csv("../data/lebensdauer.txt", "\t", 1)
    end_of_measurement = 6500
    bins = (data[0])[100:end_of_measurement] # + [data[-1] + 1]
    counts = p.ev(data[1], np.sqrt(data[1]))[100:end_of_measurement]
    p0 = [395000, 12 * 70, 1800, 170]
    res, (err, rsq) = std.odr_fit(lt_spectrum, bins, counts, p0)
    print(res)
    print((p.ev(res[1], err[1]) * time_slope).format())
    print((np.log(2) * p.ev(res[1], err[1]) * time_slope).format())
    # res, (err, rsq) = std.curve_fit(lt_spectrum, bins, counts, res)
    std.default.plt_errorbar(bins, counts, "Messdaten")
    std.default.plt_func(lt_spectrum, res, label=f"$R^2 = {round(rsq, 3)}$")
    std.default.plt_finish("Bin", "Anzahl")
