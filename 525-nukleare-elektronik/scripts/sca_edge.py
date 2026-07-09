import std
import propeller as p
import numpy as np
from scipy import special as sp


def load_data(element, side):
    fname = f"../data/{side}/{element}_fenster.txt"
    data = std.load_csv(fname, "\t", skiprows=1)
    return data[0][:], p.ev(data[1], np.sqrt(data[1]))[:]


def get_calibration(side):
    fname = f"../data/fits/{side}/energy_calibration.csv"
    fit_params = std.load_csv(fname, skiprows=1)
    return lambda x:fit_params[0] * x + fit_params[1]

def rising_edge_func(x, x0, b):
    return 0.5 * (1 + sp.erf(b * (x - x0)))


# def falling_edge_func(x, x0, b, c):
#     return 0.5 * (1 - sp.erf(b * (x - x0))) / c + c
def falling_edge_func(x, x0, b):
    return 0.5 * (1 - sp.erf(b * (x - x0)))


def window_func(x, xr, xf, br, bf, c, a, mu, sigma):
    # return rising_edge_func(x, xr, br) * (falling_edge_func(x, xf, bf) * std.area_gaussian(x, a, mu, sigma) + c)
    return rising_edge_func(x, xr, br) * falling_edge_func(x, xf, bf) * std.area_gaussian(x, a, mu, sigma) + rising_edge_func(x, xr, br) * c


p0s = {
    "na": {
        "links": [6575, 7675, 0.5, 0.5, 75, 2e6, 7300, 230],
        "rechts": [6745, 7745, 0.5, 0.1, 75, 2e6, 7300, 230]
    },
    "ba": {
        "links": None,
        "rechts": None
    }
}

def fit_edge(element, side):
    x, y = load_data(element, side)
    lb_a = 6500 if side == "links" else 6650
    ub_b = 7750 if side == "links" else 7775

    p0 = p0s[element][side]
    res, (err, rsq) = std.odr_fit(window_func, x[lb_a:ub_b], y[lb_a:ub_b], p0)

    energy_cal = get_calibration(side)
    print(p.ev(res, err))
    print(energy_cal(res[0]))
    print(energy_cal(res[1]))
    std.default.plt_errorbar(x[6500:7850], y[6500:7850])
    red_chi2 = std.reduced_chi_2(~y[lb_a:ub_b], window_func(x[lb_a:ub_b], *res), res)
    print(red_chi2)
    # red_chi2 = 1
    std.default.plt_func(window_func, res, xrange=(lb_a, ub_b), label="$\\chi^2_\\text{red}" + f"{round(red_chi2,  3)}$")
    std.default.plt_finish("Bins", "Anzahl", f"../figs/{side}/sca_{element}_window_fit.pdf")

if __name__ == "__main__":
    fit_edge("na", "links")
    fit_edge("na", "rechts")

