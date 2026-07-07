import std
import propeller as p
import numpy as np
from scipy import special as sp


def load_data(side):
    fname = f"../data/{side}/na_cfd.txt"
    data = std.load_csv(fname, "\t", skiprows=1)
    return data[0][:1500], p.ev(data[1], np.sqrt(data[1]))[:1500]



def edge_func(x, x0, a, b):
    return 0.5 * a * (1 + sp.erf(b * (x - x0)))


def fit_edge(side):
    x, y = load_data(side)
    p0 = [500, 250, 0.01]
    res, (err, rsq) = std.odr_fit(edge_func, x, y, p0)
    print(p.ev(res, err))
    std.default.plt_errorbar(x, y)
    std.default.plt_func(edge_func, res, f"$R^2 = {round(rsq, 3)}$")
    std.default.plt_finish("", "")

if __name__ == "__main__":
    fit_edge("links")
    fit_edge("rechts")
