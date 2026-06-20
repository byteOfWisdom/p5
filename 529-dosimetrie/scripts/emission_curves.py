
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


print(metadata)
plot_measurement("current_1")
plot_measurement("current_2")
std.default.plt_finish("Heizstrom $I_H$ / A", "Ionisationsstrom $I_C$ / A")


plot_measurement("voltage_1")
plot_measurement("voltage_2")
std.default.plt_finish("Beschleunigungsspannung $U_B$ / V", "Ionisationsstrom $I_C$ / A")
