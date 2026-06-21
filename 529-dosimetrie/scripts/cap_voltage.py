# %%
import std
import tomllib
import propeller as p

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["dosimetrie"]
fhandle.close()

def get_data(measurement):
    print(metadata[measurement])
    data = std.util.load_csv("../data/" + metadata[measurement]["csv_file"], skiprows=1)
    cap_voltage, amp_voltage = data[0], data[1]
    resistance = p.ev(metadata[measurement]["resistor"], 2e-2 * metadata[measurement]["resistor"])
    ion_current = amp_voltage / resistance
    return cap_voltage, amp_voltage, ion_current



def plot_measurement(measurement):
    cap_voltage, _, ion_current = get_data(measurement)
    label = f"$U_B$ = {metadata[measurement]["acceleration_voltage"] * 1e-3} kV"
    std.default.plt_errorbar(cap_voltage, ion_current, label)

# %%
# print(metadata)
for measurement in metadata.keys():
    if metadata[measurement]["element"] == "Mo":
        plot_measurement(measurement)

std.default.plt_finish("Kondensatorspannung $U_C$ / V", "Ionisationsstrom $I_C$ / A")


# %%
for measurement in metadata.keys():
    if metadata[measurement]["element"] == "Cu":
        plot_measurement(measurement)

std.default.plt_finish("Kondensatorspannung $U_C$ / V", "Ionisationsstrom $I_C$ / A")

# %%
a, b = 0, 0
for key in metadata.keys():
    u_c, u_amp, i_ion = get_data(key)
    table = {
        "$U_C$ / V": u_c,
        "$U_{amp} / V$": u_amp,
        "$I_C$ / A": i_ion
    }

    if metadata[key]["element"] == "Cu":
        a += 1
    else:
        b += 1

    fname_ending = str(a) if metadata[key]["element"] == "Cu" else str(b)

    std.print_tex_table(table, "../latex/cap_voltage" + metadata[key]["element"] + fname_ending + ".table")
