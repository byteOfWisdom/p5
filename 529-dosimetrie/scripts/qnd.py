import std
import tomllib

fhandle = open("../data/meta.toml", "rb")
metadata = tomllib.load(fhandle)["dosimetrie"]
fhandle.close()


def plot_measurement(measurement):
    print(metadata[measurement])
    data = std.util.load_csv("../data/" + metadata[measurement]["csv_file"], skiprows=1)
    print(data[0])
    print(data[1])
    std.default.plt_errorbar(data[0], data[1])


plot_measurement("messung_1")
plot_measurement("messung_2")
plot_measurement("messung_3")
plot_measurement("messung_4")
std.default.plt_finish("x", "y")
