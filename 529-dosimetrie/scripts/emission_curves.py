# %%
import std
import tomllib
import propeller as p
from matplotlib import pyplot as plt


fhandle = open("../data/meta.toml", "rb")
all = tomllib.load(fhandle)
metadata = all["emission"]
weather = all["roominfo"]
fhandle.close()

air_volume = 125e-6 # cubic meter?
rho_0 = 1293
pressure = p.ev(weather["air_pressure"], 3)
t0 = 273
t = p.ev(weather["air_temp"], 2.0)
p0 = 1013
air_density = rho_0 * (t0 / t) * (pressure / p0)
air_mass = air_volume * air_density


def plot_measurement(measurement):
    print(metadata[measurement])
    data = std.util.load_csv("../data/" + metadata[measurement]["csv_file"], skiprows=1)
    x_var, amp_voltage = data[0], data[1]
    resistance = p.ev(metadata[measurement]["resistor"], 2e-2 * metadata[measurement]["resistor"])
    ion_current = amp_voltage / resistance
    ion_dose_rate = ion_current / air_mass
    label = metadata[measurement]["element"]
    std.default.plt_errorbar(x_var, ion_dose_rate, label)
    return x_var, ion_current, ion_dose_rate


# %%
heating_current_cu, ion_current_cu, ion_dose_rate_cu = plot_measurement("current_1")
cu_res, (cu_err, cu_rsq) = std.curve_fit(std.linear, heating_current_cu, ion_dose_rate_cu)
std.default.plt_func(std.linear, cu_res, f"$R^2= {round(cu_rsq, 3)}$")
heating_current_mo, ion_current_mo, ion_dose_rate_mo = plot_measurement("current_2")
mo_res, (mo_err, mo_rsq) = std.curve_fit(std.linear, heating_current_mo, ion_dose_rate_mo)
std.default.plt_func(std.linear, mo_res, f"$R^2 = {round(mo_rsq, 3)}$")

dose_y = plt.gca().secondary_yaxis("right",functions=(lambda x: 32.4 * x, lambda x: x / 32.4))
dose_y.set_ylabel("Äquivalentsdosisleistung / $Svs^{-1}$")
std.default.plt_finish("Heizstrom $I_H$ / A", "mittlere Ionisationsdosisleistung <j> / $Akg^{-1}$")

cu_current_table = {
    "I_H / mA": heating_current_cu,
    "I_C / A": ion_current_cu,
    "<j> / $Akg^{-1}$": ion_dose_rate_cu,
}
std.print_tex_table(cu_current_table, "../latex/cu_current.table")

mo_current_table = {
    "I_H / mA": heating_current_mo,
    "I_C / A": ion_current_mo,
    "<j> / $Akg^{-1}$": ion_dose_rate_mo,
}
std.print_tex_table(mo_current_table, "../latex/mo_current.table")

# %%
a, b, c = plot_measurement("voltage_1")
voltage_table = {
    "U_H / V": a,
    "I_C / A": b,
    "<j> / $Akg^{-1}$": c,
}

std.print_tex_table(voltage_table, "../latex/cu_voltage.table")

a, b, c = plot_measurement("voltage_2")
voltage_table = {
    "U_H / V": a,
    "I_C / A": b,
    "<j> / $Akg^{-1}$": c,
}

std.print_tex_table(voltage_table, "../latex/mo_voltage.table")


dose_y = plt.gca().secondary_yaxis("right",functions=(lambda x: 32.4 * x, lambda x: x / 32.4))
dose_y.set_ylabel("Äquivalentsdosisleistung / $Svs^{-1}$")
std.default.plt_finish("Beschleunigungsspannung $U_B$ / V",  "mittlere Ionisationsdosisleistung <j> / $Akg^{-1}$")

