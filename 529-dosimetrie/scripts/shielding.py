import std
import propeller as p
import tomllib
import numpy as np

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["shielding"]
fhandle.close()

thickness_lut = np.array([0, 0.5, 1., 1.5, 2., 2.5, 3.])

# a = p.from_string("0.4882(35)")
# b = p.from_string("0.01080(22)")


def current_corrected_countrate(current, countrate):
    a = p.ev(0.4882, 0.0035)
    b = p.ev(0.01080, 0.00022)
    i_max = max(current)
    lin_corr = [(a * i_max + b) / (a * i + b) for i in current]
    # return countrate * lin_corr
    current = current / max(current)
    return countrate / current # TODO: is this good enough?

ccc = current_corrected_countrate

def absorber(x, u):
    return np.exp(- x * u)


def process_measurement(params):
    print(params)
    data = std.load_csv("../data/" + params["csv_file"], skiprows=1)
    angle, current, countrate = data[0], data[1], data[2]
    key = np.argsort(angle)
    current = current[key]
    current = p.ev(current, 0.01)
    countrate = countrate[key]
    countrate = p.ev(countrate, np.sqrt(countrate))
    actual_rate = ccc(current, countrate)
    # actual_rate = countrate / current

    thickness = thickness_lut
    transmission = actual_rate / actual_rate[0]

    res, (err, rsq) = std.curve_fit(absorber, thickness[1:], transmission[1:])
    std.default.plt_errorbar(thickness, transmission, "mit Filter" if params["filter"] else "ohne Filter")
    std.default.plt_func(absorber, res, f"$R^2 = {rsq}$", (0, 3))

    return p.ev(res, err)



res_a = process_measurement(metadata["messung_1"])
res_b = process_measurement(metadata["messung_2"])
print(res_a)
print(res_b)
std.default.plt_func(absorber, [0.929], "refrence", (0, 3))
std.default.plt.yscale("log")
std.default.plt_finish("Dicke / mm", "Transmissivität / 1")
