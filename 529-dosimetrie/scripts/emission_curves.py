
import std
import tomllib
import propeller as p

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["emission"]
fhandle.close()


def plot_measurement(measurement):
    print(metadata[measurement])
    data = std.util.load_csv("../data/" + metadata[measurement]["csv_file"], skiprows=1)
    x_var, amp_voltage = data[0], data[1]
    resistance = p.ev(metadata[measurement]["resistor"], 2e-2 * metadata[measurement]["resistor"])
    ion_current = amp_voltage / resistance
    label = metadata[measurement]["element"] 
    std.default.plt_errorbar(x_var, ion_current, label)
    return x_var, ion_current


print(metadata)
heating_current_cu, ion_current_cu = plot_measurement("current_1")
cu_res, (cu_err, cu_rsq) = std.curve_fit(std.linear, heating_current_cu, ion_current_cu)
std.default.plt_func(std.linear, cu_res, f"$R^2= {round(cu_rsq, 3)}$")
heating_current_mo, ion_current_mo = plot_measurement("current_2")
mo_res, (mo_err, mo_rsq) = std.curve_fit(std.linear, heating_current_mo, ion_current_mo)
std.default.plt_func(std.linear, mo_res, f"$R^2 = {round(mo_rsq, 3)}$")
std.default.plt_finish("Heizstrom $I_H$ / A", "Ionisationsstrom $I_C$ / A")


plot_measurement("voltage_1")
plot_measurement("voltage_2")
std.default.plt_finish("Beschleunigungsspannung $U_B$ / V", "Ionisationsstrom $I_C$ / A")
