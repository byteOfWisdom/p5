import numpy as np
import std
from std import odr_fit
import propeller as p
from matplotlib import pyplot as plt
from sys import argv


def get_data(file):
    volts, current_volt = std.load_csv(file, delimiter = ",", skiprows = 1)
    volts = volts*1e3
    volts = p.ev(volts, volts*1e-2) #in Volt
    current_volt = current_volt*1e-3 #in Volt
    return volts, current_volt

def plot_data(file, color="tab:green"):
    volts, current_volt = get_data(file)
    resist = 1e6
    resist = p.ev(resist, resist*1e-2)
    currents = current_volt / resist
    #print(currents)
    _, v_err = p.ve(volts)
    _, c_err = p.ve(currents)
    #label = str(file.removeprefix("../data/").removesuffix(".csv").replace("_"," "))

    if "ohne" in str(file):
        label = "Messung ohne Präperat"
        color = "red"
    else:
        label = "Messung mit Präperat"
    plt.errorbar(~volts, ~currents, yerr = c_err, xerr = v_err, label = label, color = color, **std.default.error_bar_def)
    return

def fit_exp(file):
    volts, curr_volts = get_data(file)
    resist = 1e6
    resist = p.ev(resist, resist*3e-2)
    #_, r_err = p.ve(resist)
   #print(r_err)
    currents = curr_volts / resist
    params, (std, goodness) = odr_fit(exponential, volts, currents)
    print("fitfunktion: c*np.exp(a*x)")
    #print("a, c")
    for i in range(len(params)):
        print("a: ") if i == 0 else print("c: ")
        fit = p.ev(params[i], std[i])
        print(fit.format())
    x_ax = np.linspace(min(volts),max(volts),1000)
    plt.plot(x_ax, exponential(x_ax, *params), linestyle = "-", linewidth = 0.8, label = f"Anpassungsfunktion, $R^2$={str(round(goodness,3)).replace(".",",")}")
    return


def exponential(x,a,c):
     return c*np.exp(a*(x))


def plot(file1, file2, saved=False):
    plot_data(file1)
    plot_data(file2)
    fit_exp(file1)
    plt.legend(loc="best")
    std.default.plt_pretty("Spannung / V", "Strom / A")
    if saved:
        plt.savefig("../figs/strom_vergleich.pdf")
        plt.clf()
    else:
        plt.show()
    # #only without
    plot_data(file2)
    plt.legend(loc="best")
    std.default.plt_pretty("Spannung / V", "Strom / A")
    if saved:
        plt.savefig("../figs/strom_ohne.pdf")
    else:
        plt.show()

    #plot_data(file1)
    #fit_exp(file1)
    # plt.legend(loc="best")
    # std.default.plt_pretty("Spannung / V", "Strom/ A")
    # if saved:
    #     plt.savefig("../figs/strom_fit.pdf")
    # else:
    #     plt.show()
    return

def main():
    plot(argv[1], argv[2], saved=False)
    return


if __name__ == "__main__":
    main()
