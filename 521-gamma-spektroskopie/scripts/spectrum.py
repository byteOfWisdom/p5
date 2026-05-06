from sys import argv
import numpy as np
import spectrum_fit
import std
from matplotlib import pyplot as plt


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    ug_data = np.transpose(np.loadtxt(argv[2]))

    plt.plot(data[0], data[1], linewidth=0.5)
    out_name = argv[1][:-3].split("/")[-1] + "pdf"
    plt.title(out_name[:-4])
    std.default.plt_finish("Kanal", "Count", "figs/" + out_name)
    return None

    channels = data[0]
    hits = data[1] - ug_data[1]
    lines = spectrum_fit.analyze_spectrum(channels[:-1], hits[:-1], argv[3] if len(argv) > 3 else False)
    if len(argv) > 3:
        std.print_tex_table(lines, argv[3] + ".table")
        std.print_csv_table(lines, argv[3] + ".csv")
    else:
        print(lines)
    return None


if __name__ == "__main__":
    main()
