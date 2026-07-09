import std
import numpy as np
import propeller as p
import time_cal
import scipy.special as sp
import iminuit as im


def lt_spectrum(t, amplitude, tau, t0,  sigma):
    dt = t - t0
    # b = ((dt / sigma) - (sigma / tau)) / np.sqrt(2)
    b = (tau * dt - sigma**2) / (np.sqrt(2) * sigma * tau)
    # return b
    # d = (sigma ** 2 + tau * t0) / (np.sqrt(2) * sigma * tau)
    c = amplitude * (1 + sp.erf(b)) / (2 * tau)
    return c * np.exp(((sigma**2) - (2 * tau * dt)) / (2 * tau ** 2))


def iminuit_fit(bins, counts, p0):
    # least_squares = im.cost.LeastSquares(bins, ~counts, p.ve(counts)[1], lt_spectrum)
    # m = im.Minuit(least_squares, *p0)
    binned_nll = im.cost.ExtendedBinnedNLL(~counts, np.append(bins, bins[-1] + 1), lt_spectrum)
    m = im.Minuit(binned_nll, *p0)
    m.simplex()
    m.migrad()
    m.hesse()
    print(m)
    return m.values


if __name__ == "__main__":
    time_cal_params = time_cal.invert(*time_cal.time_cal())
    time_slope = time_cal_params[0]
    data = std.load_csv("../data/lebensdauer.txt", "\t", 1)
    end_of_measurement = 6500
    bins = (data[0])[100:end_of_measurement] # + [data[-1] + 1]
    counts = p.ev(data[1], np.sqrt(data[1]))[100:end_of_measurement]
    p0 = [395000, 12 * 70, 1800, 170]
    # res = iminuit_fit(bins, counts, p0)
    res, (err, rsq) = std.odr_fit(lt_spectrum, bins, counts, p0)
    print(res)
    tau = p.ev(res[1], err[1]) * time_slope
    hlt = np.log(2) * tau
    print(tau.format())
    print(hlt.format())
    result_str = f"{std.si_string("\\tau", tau, "\\nano\\second")} und {std.si_string("T_\\frac{1}{2}", hlt, "\\nano\\second")}"
    std.write_file("../data/fits/lifetime.txt", result_str)

    params = p.ev(res, err)
    fit_res_tbl = {
        "$I$ / Anzahl": [params[0]],
        "$\\tau$ / Bins": [params[1]],
        "$t_0$ / Bins": [params[2]],
        "$\\sigma$ / Bins": [params[3]],
    }
    std.print_tex_table(fit_res_tbl, "../latex/lifetime_fit.table")

    chi2 = std.reduced_chi_2(~counts, lt_spectrum(bins, *res), res)
    # res, (err, rsq) = std.curve_fit(lt_spectrum, bins, counts, res)
    std.default.plt_errorbar(bins, counts, "Messdaten")
    std.default.plt_func(lt_spectrum, res, label="$\\chi^2_{red} = " + f"{round(chi2, 3)}$")
    std.default.plt_finish("Bin", "Anzahl", "../figs/lifetime_plot.pdf")
