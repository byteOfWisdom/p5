import numpy as np
from matplotlib import pyplot as plt
import std
from sys import argv
import propeller as p
from scipy.optimize import curve_fit


def get_data(file):
     # get lit vals
    literature = lit_values(str(file))
    #print("lit vals:", literature)
    ref_vals = np.array(literature)

    #get measurement data
    try:
        measurement = std.util.load_csv("../figs/"+file+"_nai_bin.csv", skiprows=1)
    except FileNotFoundError:
        return print("invalid file name!")
    mu_bins = measurement[1]
    #print(mu_bins)

    # select correct lines to work with
    if file == "eu":
        #only use certain elements
        mu_bins = np.append(mu_bins[4:7],mu_bins[8:])
    #print("after slicing: ",mu_bins)
    peak_lines = np.array([mu_bins,ref_vals])
    # print(peak_lines)
    #print("measured mu:",peak_lines[0])
    #print("lit vals:",peak_lines[1])
    return peak_lines[0], peak_lines[1]


def lit_values(element):
    cobalt = [1173.2,1332.5]
    caesium = [661.7, 283.5]
    europium = [121.7817, 244.6974, 344.2785, 778.9045, 964.057, 1085.837, 1408.013]
    if str(element) == "co":
        return cobalt
    if str(element) == "cs":
        return caesium
    if str(element) == "eu":
        return europium
    return "error with the lit vals!"


def plot_data(save_to=False):
    x_meas_eu, y_lit_eu = get_data("eu")
    x_meas_co, y_lit_co = get_data("co")
    x_meas_cs, y_lit_cs = get_data("cs")

    x_temp = np.append(x_meas_co,x_meas_cs[:-1])
    x_meas = np.append(x_temp, x_meas_eu)
    x_meas_vals, x_meas_errs = p.ve(x_meas)
    #print(x_meas_errs)
    y_temp = np.append(y_lit_co, y_lit_cs[:-1])
    y_lit = np.append(y_temp, y_lit_eu)
    #print(type(y_lit[0]))

    #print(len(x_meas), len(y_lit))
    fit, cov = curve_fit(lambda x, a, b: a*x+b, x_meas_vals, y_lit*1e3, bounds=([0,0],[1000,1e10]))
    err = np.sqrt(np.diag(cov))

    params = []
    for i in range(len(fit)):
        params.append(p.ev(fit[i],err[i]))

    #print(y_lit*1e3)
    #print(std.linear(x_meas,*fit))
    goodness = std.math_stuff.goodness_of_fit(y_lit*1e3, std.linear(x_meas_vals,*fit))
    #print(fit)
    #print(err)
    x_ax = np.linspace(0, max(~x_meas), 500)
    #print(fit[0]*90+fit[1])

    plt.errorbar(x_meas_vals, y_lit*1e3, xerr=x_meas_errs,label="zugeordnete Linien (NaI-Detektor)",color="tab:red", **std.default.error_bar_def)
    #plt.scatter(x_meas_vals, y_lit*1e3, color="tab:blue", marker="x", linewidths=0.7)

    plt.plot(x_ax,std.linear(x_ax,*fit), label=f"Anpassungsgerade, $R^2$={round(goodness,3)}",color="tab:blue")

    plt.legend(loc="best", fontsize="small")
    std.default.plt_pretty("Kanalnr. / [Anzahl]","Energie / eV")

    print(params[0].format(),params[1].format())

    if save_to:
        plt.savefig("../figs/nai_calibration.pdf")

    else:
        plt.show()
    return fit, err



def main():
    #get_data("co")
    fit, err = plot_data(save_to=True)
    #print(fit, err)
    #print(goodness)
    return

if __name__ == "__main__":
    main()
