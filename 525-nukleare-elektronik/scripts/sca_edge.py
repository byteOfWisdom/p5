import std
import propeller as p
import numpy as np
from scipy import special as sp


def load_data(element, side):
    fname = f"../data/{side}/{element}_fenster.txt"
    if element == "ba":
        fname = "../data/links/start_sca_356.txt" if side == "links" else "../data/rechts/stop_sca_81.txt"
    data = std.load_csv(fname, "\t", skiprows=1)
    print(data)
    return data[0][:], p.ev(data[1], np.sqrt(data[1]))[:]


def get_calibration(side):
    fname = f"../data/fits/{side}/energy_calibration.csv"
    fit_params = std.load_csv(fname, skiprows=1)
    return lambda x:fit_params[0][0] * x + fit_params[1][0]

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
        "links": [4500, 5270, 0.5, 0.5, 5, 2e6, 4951, 230],
        "rechts": [1120, 1340, 0.5, 0.5, 5, 3e6, 1235, 100]
    }
}

bounds = {
    "na": {
        "links": (6500, 7850),
        "rechts": (6500, 7850)
    },
    "ba": {
        "links": (4250, 5500),
        "rechts": (1000, 1500)
    }
}

ublb = {
    "na": {
        "links": (6500, 7750),
        "rechts": (6650, 7775)
    },
    "ba": {
        "links": (4250, 5500),
        "rechts": (1000, 1400)
    }
}

def only_plot(element, side):
    x, y = load_data(element, side)
    std.default.plt_errorbar(x, y)
    std.default.plt_finish("Bins", "Anzahl")



def fit_edge(element, side):
    x, y = load_data(element, side)
    # lb_a = 6500 if side == "links" else 6650
    # ub_b = 7750 if side == "links" else 7775
    lb_a, ub_b = ublb[element][side]

    p0 = p0s[element][side]
    res, (err, rsq) = std.odr_fit(window_func, x[lb_a:ub_b], y[lb_a:ub_b], p0)

    param_names = ["$x_0$ / Bin", "$x_1$ / Bin", "$a_0$ / $\\text{Bin}^{-1}$", "$a_1$ / $\\text{Bin}^{-1}$", "c / Bin", "A / Anzahl", "$\\mu$ / Bin", "$\\sigma$ / Bin"]

    energy_cal = get_calibration(side)
    print(p.ev(res, err))
    params = p.ev(res, err)
    table = {
        "Parameter": param_names,
        "Wert": params
    }
    # ---------
    std.print_tex_table(table, f"../latex/{side}_{element}_sca_window.table")

    print(energy_cal(params[0]))
    print(energy_cal(params[1]))
    lower_edge = energy_cal(params[0])
    upper_edge = energy_cal(params[1])
    value_string = f"{std.si_string("x_0", lower_edge, "\\kilo\\eV")} und {std.si_string("x_1", upper_edge, "\\kilo\\eV")}"

    # ---------
    std.write_file(f"../data/fits/{side}/{element}_sca_edges.txt", value_string)
    std.write_file(f"../data/fits/{side}/{element}_sca_interval.txt", std.si_string("\\Delta E", upper_edge - lower_edge, "\\kilo\\eV"))

    plot_from, plot_to = bounds[element][side]
    std.default.plt_errorbar(x[plot_from:plot_to], y[plot_from:plot_to])
    red_chi2 = std.reduced_chi_2(~y[lb_a:ub_b], window_func(x[lb_a:ub_b], *res), res)
    print(red_chi2)
    # red_chi2 = 1
    std.default.plt_func(window_func, res, xrange=(lb_a, ub_b), label="$\\chi^2_\\text{red}" + f"{round(red_chi2,  3)}$")

    # ----------
    std.default.plt_finish("Bins", "Anzahl", f"../figs/{side}/sca_{element}_window_fit.pdf")
    # std.default.plt_finish("Bins", "Anzahl")

if __name__ == "__main__":
    fit_edge("na", "links")
    fit_edge("na", "rechts")
    fit_edge("ba", "links")
    fit_edge("ba", "rechts")

