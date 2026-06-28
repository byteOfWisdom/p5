import std
import propeller as p
import tomllib
import numpy as np
import xraydb

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["shielding"]
fhandle.close()

thickness_lut = np.array([0, 0.5, 1., 1.5, 2., 2.5, 3.])

# a = p.from_string("0.4882(35)")
# b = p.from_string("0.01080(22)")


def current_corrected_countrate(current, countrate):
    # a = p.ev(0.4882, 0.0035)
    # b = p.ev(0.01080, 0.00022)
    # lin_corr = [(a * i_max + b) / (a * i + b) for i in current]
    # return countrate * lin_corr
    # current = current / max(current)
    return countrate / current # TODO: is this good enough?


ccc = current_corrected_countrate


def absorber(x, u):
    return np.exp(- x * u)


table_file = None

def load_and_normalize(file):
    data = std.load_csv("../data/" + file, skiprows=1)
    angle, current, countrate = data[0], data[1], data[2]

    # sort in ascending order
    key = np.argsort(angle)
    current = current[key]
    countrate = countrate[key]

    current = p.ev(current, 0.01)
    countrate = p.ev(countrate, np.sqrt(countrate))
    actual_rate = ccc(current, countrate)

    transmission = actual_rate / actual_rate[0]
    table = {
        "I_H / mA": current,
        "Zählrate / $s^{-1}$": countrate,
        "korrigierte Zählrate / $s^{-1}$": actual_rate,
        "Transmissivität": transmission
    }

    std.print_tex_table(table, table_file)
    return transmission


def process_measurement(params):
    print(params)
    transmission = load_and_normalize(params["csv_file"])

    thickness = p.ev(thickness_lut, 0.05)

    res, (err, rsq) = std.curve_fit(absorber, thickness[1:], transmission[1:])
    filter_str = "mit Filter" if params["filter"] else "ohne Filter"
    std.default.plt_errorbar(thickness, transmission, filter_str)
    std.default.plt_func(absorber, res, f"{filter_str} $R^2 = {round(rsq, 3)}$", (0, 3))

    return p.ev(res, err)


if __name__ == "__main__":
    table_file = "../latex/abschirmung_ohne.table"
    res_a = process_measurement(metadata["messung_1"])
    print(res_a)
    table_file = "../latex/abschirmung_mit.table"
    res_b = process_measurement(metadata["messung_2"])
    print(res_b)
    literature_mu = 0.1 * xraydb.material_mu("Al", 20e3)
    print(literature_mu)
    std.default.plt_func(absorber, [literature_mu], "refrence", (0, 3))
    std.default.plt.yscale("log")
    std.default.plt_finish("Dicke / mm", "Transmissivität / 1")
