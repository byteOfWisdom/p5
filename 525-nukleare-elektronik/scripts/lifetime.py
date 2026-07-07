import std
import numpy as np
import propeller as p
import time_cal
import scipy.special as sp
import iminuit as im


def lt_spectrum(t, amplitude, tau, t0,  sigma):
    dt = t - t0
    b = ((dt / sigma) - (sigma / tau)) / np.sqrt(2)
    c = amplitude * (1 + sp.erf(b)) / (2 * tau)
    return c * np.exp(((sigma**2) - (2 * tau * dt)) / (2 * tau ** 2))


if __name__ == "__main__":
    data = std.load_csv("../data/lebensdauer.txt", "\t", 1)
    bins = data[0] # + [data[-1] + 1]
    counts = p.ev(data[1], np.sqrt(data[1]))
    p0 = [~max(counts), 1, 2000, 100]
    res, (err, rsq) = std.odr_fit(lt_spectrum, bins, counts, p0)
    print(res)
    # c = im.cost.BinnedNLL(data[1], bins, lt_spectrum)
    # m = im.Minuit(c)
    # m.migrad()
    print(lt_spectrum(bins, *res))
    std.default.plt_errorbar(bins, counts)
    std.default.plt_func(lt_spectrum, res)
    std.default.plt_finish("Bin", "Anzahl")
