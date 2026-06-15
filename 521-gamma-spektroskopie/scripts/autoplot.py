import std
import numpy as np
from matplotlib import pyplot as plt

files = [
    "co_ge.txt",
    "co_nai.txt",
    "cs_ge.txt",
    "cs_nai_fixed.txt",
    "eu_ge.txt",
    "eu_nai.txt",
    "langzeit_mitprobe_ge.txt",
    "langzeit_ohneprobe_ge.txt",
    "undergrd_ge.txt",
    "undergrd_nai.txt",
]

for file in files:
    print(file)
    data = std.load_csv("data/" + file)
    x = data[0][1:-1]
    y = data[1][1:-1]
    yerr = np.sqrt(y)

    plt.cla()
    plt.errorbar(x, y, yerr, label="Messdaten", **std.default.error_bar_def)
    std.default.plt_finish("Kanal", "Anzahl", save_to= "figs/" + file[:-4] + "_with_eb.pdf")
