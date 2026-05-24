import numpy as np
from matplotlib import pyplot as plt
#import scipy
import std
#import propeller as p
from sys import argv

def get_data(file):
    times, currents = np.genfromtxt(file, delimiter=",",skip_header=18, usecols=(3,4),autostrip=True, unpack=True)
    #print(times[:10])
    #print(currents[:10])
    #currents.strip()
    return times, currents

def plot_data(file):
    times, currents = get_data(file)
    #print(times[:10])
    #print(currents[:10])
    t_err = []
    c_err = []
    for t in range(len(times)):
        t_err.append(2*1e-9) #set times err to time div/2
        c_err.append(abs(currents[t]*3e-2)) #set current errs to 3%
    plt.errorbar(times, currents, xerr=t_err, yerr=c_err, label="Messdaten", **std.default.error_bar_def)
    std.default.plt_pretty("Zeit / s", "Spannung / V")
    plt.legend()
    plt.show()
    return


def main():
    plot_data(argv[1])
    return

if __name__ == "__main__":
    main()
