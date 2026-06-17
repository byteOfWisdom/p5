import std
from matplotlib import pyplot as plt
import tomllib
from sys import argv

fhandle = open("../data/plate_cap.toml", "rb")
metadata = tomllib.load(fhandle)
fhandle.close()



data = std.util.load_csv("../data/" + metadata["dosimetrie"][argv[1]]["csv_file"], skiprows=1)
plt.scatter(data[0], data[1])
plt.show()
