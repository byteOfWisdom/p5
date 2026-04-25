from sys import argv
from matplotlib import pyplot as plt
import numpy as np


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    plt.plot(data[0], data[1])
    plt.show()


if __name__ == "__main__":
    main()
