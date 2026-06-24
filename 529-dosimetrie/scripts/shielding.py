import std
import propeller as p
import tomllib
import numpy as np

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["shielding"]
fhandle.close()

thickness_lut = np.array([0, 0.5, 1., 1.5, 2., 2.5, 3.])

def current_corrected_countrate(current, countrate):
    return countrate / current # TODO: is this good enough?

ccc = current_corrected_countrate

def absorber(x, u):
    return np.exp(- x * u)

def process_measurement(params):
    print(params)
    data = std.load_csv("../data/" + params["csv_file"], skiprows=1)
    angle, current, countrate = data[0], data[1], data[2]
    countrate = p.ev(countrate, np.sqrt(countrate))
    actual_rate = ccc(current, countrate)

    thickness = thickness_lut[np.argsort(angle)]
    transmission = actual_rate / max(actual_rate)

    res, (err, rsq) = std.curve_fit(absorber, thickness, transmission)
    std.default.plt_errorbar(thickness, transmission, "mit Filter" if params["filter"] else "ohne Filter")
    std.default.plt_func(absorber, res, f"$R^2 = {rsq}$", (0, 3))

    return p.ev(res, err)



process_measurement(metadata["messung_1"])
process_measurement(metadata["messung_2"])
std.default.plt.yscale("log")
std.default.plt_finish("Dicke / mm", "Transmissivität / 1")
