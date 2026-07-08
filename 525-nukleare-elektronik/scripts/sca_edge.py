import std
import propeller as p
import numpy as np
from scipy import special as sp


def load_data(side):
    fname = f"../data/{side}/na_fenster.txt"
    data = std.load_csv(fname, "\t", skiprows=1)
    return data[0][:], p.ev(data[1], np.sqrt(data[1]))[:]



def rising_edge_func(x, x0, b):
    return 0.5 * (1 + sp.erf(b * (x - x0)))


# def falling_edge_func(x, x0, b, c):
#     return 0.5 * (1 - sp.erf(b * (x - x0))) / c + c
def falling_edge_func(x, x0, b):
    return 0.5 * (1 - sp.erf(b * (x - x0)))


def window_func(x, xr, xf, br, bf, c, a, mu, sigma):
    # return rising_edge_func(x, xr, br) * (falling_edge_func(x, xf, bf) * std.area_gaussian(x, a, mu, sigma) + c)
    return rising_edge_func(x, xr, br) * falling_edge_func(x, xf, bf) * std.area_gaussian(x, a, mu, sigma) + rising_edge_func(x, xr, br) * c


def fit_edge(side):
    x, y = load_data(side)
    lb_a = 6500 if side == "links" else 6650
    ub_b = 7750 if side == "links" else 7775

    p0 = [6575, 7675, 0.5, 0.5, 75, 2e6, 7300, 230]
    if side != "links":
        p0 = [6745, 7745, 0.5, 0.1, 75, 2e6, 7300, 230]
    res, (err, rsq) = std.odr_fit(window_func, x[lb_a:ub_b], y[lb_a:ub_b], p0)

    print(p.ev(res, err))
    std.default.plt_errorbar(x[6500:7850], y[6500:7850])
    red_chi2 = std.reduced_chi_2(~y[lb_a:ub_b], window_func(x[lb_a:ub_b], *res), res)
    print(red_chi2)
    # red_chi2 = 1
    std.default.plt_func(window_func, res, xrange=(lb_a, ub_b), label=round(red_chi2,  3))
    std.default.plt_finish("Bins", "Anzahl")

if __name__ == "__main__":
    fit_edge("links")
    fit_edge("rechts")

