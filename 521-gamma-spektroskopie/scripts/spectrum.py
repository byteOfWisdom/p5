from sys import argv
import numpy as np
import spectrum_fit
import std


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    ug_data = np.transpose(np.loadtxt(argv[2]))

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
