
import std
import tomllib
import propeller as p

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["dosimetrie"]
fhandle.close()


def plot_measurement(measurement):
    print(metadata[measurement])
    data = std.util.load_csv("../data/" + metadata[measurement]["csv_file"], skiprows=1)
    cap_voltage, amp_voltage = data[0], data[1]
    resistance = p.ev(metadata[measurement]["resistor"], 2e-2 * metadata[measurement]["resistor"])
    ion_current = amp_voltage / resistance
    label = f"$U_B$ = {metadata[measurement]["acceleration_voltage"] * 1e-3} kV" 
    std.default.plt_errorbar(cap_voltage, ion_current, label)


print(metadata)
# std.default.plt_finish("Kondensatorspannung $U_C$ / V", "Ionisationsstrom $I_C$ / A")
