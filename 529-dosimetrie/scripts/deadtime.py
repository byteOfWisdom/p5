import std
import propeller as p
import tomllib
import numpy as np

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["totzeit"]
fhandle.close()


def non_paralysing(x, a, tau):
    m = tau + (1 / (a * x))
    return 1 / m
    # return a * x / (1 - a * tau * x)


# def non_paralysing(x, a, tau):
#     m = tau + (1 / (a * x))
#     return 1 / m
#     # return a * x / (1 - a * tau * x)


def hybrid(x, a, tau_p, tau_np):
    lcc = std.linear(x, a, 0)
    # exp_term = a * x * np.exp(- x * a * tau_p)
    exp_term = lcc * np.exp(- lcc * tau_p)
    denom_term = 1 + lcc * tau_np
    return exp_term / denom_term


def paralysing(x, a, tau):
    lcc = std.linear(x, a, 0)
    # return a * x * np.exp(- x * a * tau)
    return lcc * np.exp(- lcc * tau)


print(metadata)

data = std.util.load_csv("../data/" + metadata["csv_file"], skiprows=1)
emission_current = data[0]
countrate = p.ev(data[1], np.sqrt(data[1]))
std.default.plt_errorbar(emission_current, countrate, "gemessene Zählraten")
std.default.plt_finish("Heizstrom $I_C$ / mA", "Ionendosisleistung / $Akg^{-1}$")

table = {
    "I_E / mA": emission_current,
    "Zählrate / $s^{-1}$": countrate,
}

# std.print_tex_table(table, "../latex/gmt_deadtime.table")

res_np, (err_np, rsq_np) = std.curve_fit(non_paralysing, emission_current, countrate, p0=[6e5, 6e-5])
res_p, (err_p, rsq_p) = std.curve_fit(paralysing, emission_current, countrate)
res_h, (err_h, rsq_h) = std.curve_fit(hybrid, emission_current, countrate, p0=[res_p[0], res_p[-1], res_np[-1]])

print(res_np)
print(res_p)
print(res_h)

std.default.plt_errorbar(emission_current, countrate)
std.default.plt_func(non_paralysing, res_np, label=f"nicht paraylsierend, $R^2 = {rsq_np}$", xrange=(0.01, None))
std.default.plt_func(paralysing, res_p,label=f"paraylsierend, $R^2 = {rsq_p}$", xrange=(0.01, None))
std.default.plt_func(hybrid, res_h, label=f"hybrid, $R^2 = {rsq_h}$", xrange=(0.01, None))
std.default.plt_finish("Heizstrom $I_C$ / mA", "gemessene Zählrate / $s^{-1}$")


countrate = countrate[emission_current <= 0.2]
emission_current = emission_current[emission_current <= 0.2]
# dose_rate = current_to_countrate(emission_current)#[emission_current <= 0.2]



res_np, (err_np, rsq_np) = std.curve_fit(non_paralysing, emission_current, countrate, p0=[6e5, 6e-5])
res_p, (err_p, rsq_p) = std.curve_fit(paralysing, emission_current, countrate)
res_h, (err_h, rsq_h) = std.curve_fit(hybrid, emission_current, countrate, p0=[res_p[0], res_p[-1], res_np[-1]])

print(res_np)
print(res_p)
print(res_h)

std.default.plt_errorbar(emission_current, countrate)
std.default.plt_func(non_paralysing, res_np, label=f"nicht paraylsierend, $R^2 = {rsq_np}$", xrange=(0.01, None))
std.default.plt_func(paralysing, res_p,label=f"paraylsierend, $R^2 = {rsq_p}$", xrange=(0.01, None))
std.default.plt_func(hybrid, res_h, label=f"hybrid, $R^2 = {rsq_h}$", xrange=(0.01, None))
std.default.plt_finish("Heizstrom $I_C$ / mA", "gemessene Zählrate / $s^{-1}$")
