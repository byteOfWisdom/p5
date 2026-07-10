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


def get_calibration(side):
    fname = f"../data/fits/{side}/energy_calibration.csv"
    fit_params = std.load_csv(fname, skiprows=1)
    return lambda x:fit_params[0][0] * x + fit_params[1][0]


def fit_edge(side):
    x, y = load_data(side)
    p0 = [500, 250, 0.01]
    res, (err, rsq) = std.odr_fit(edge_func, x, y, p0)
    print(p.ev(res, err))
    params = p.ev(res, err)
    energy_cal = get_calibration(side)
    fit_param_string = f"{std.si_string("x_0", params[0], "\\bin")}, {std.si_string("x_0", energy_cal(params[0]), "\\kilo\\eV")}, {std.si_string("a", params[1], "\\count")} und {std.si_string("b", params[2], "\\unity")}"
    std.write_file(f"../data/fits/{side}/cfd_edge_fit.txt", fit_param_string)
    std.default.plt_errorbar(x, y)
    std.default.plt_func(edge_func, res, f"$R^2 = {round(rsq, 3)}$")
    std.default.plt_finish("Bin", "Anzahl", f"../figs/{side}/cfd_edge.pdf")


if __name__ == "__main__":
    fit_edge("links")
    fit_edge("rechts")
