from sys import argv
from matplotlib import pyplot as plt
import numpy as np
import spectrum_fit
import std


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    lines, goodness = spectrum_fit.decomp_spectrum(data[0][:-1], data[1][:-1], lambda x, a, b: a * x + b, 2, argv[2] if len(argv) > 2 else False)
    print(lines)
    if len(argv) > 3:
        std.print_tex_table(lines, argv[3])
    return

if __name__ == "__main__":
    main()
