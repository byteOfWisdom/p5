from sys import argv
from matplotlib import pyplot as plt
import numpy as np
import spectrum_fit

def main():
    data = np.transpose(np.loadtxt(argv[1]))
    lines = spectrum_fit.decomp_spectrum(data[0], data[1], lambda x, a, b: a * x + b)
    print(lines)

if __name__ == "__main__":
    main()
