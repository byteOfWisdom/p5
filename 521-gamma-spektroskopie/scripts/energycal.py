import numpy as np
from matplotlib import pyplot as plt
import std
from sys import argv
import propeller as p
from scipy.optimize import curve_fit


def get_data(file,detector):
    if detector == "nai":
        # get lit vals
        literature = lit_values_nai(str(file))
        #print("lit vals:", literature)
        ref_vals = np.array(literature)

        #get measurement data
        try:
            measurement = std.util.load_csv("../figs/"+file+"_nai_bin.csv", skiprows=1)
        except FileNotFoundError:
            return print("invalid file name!")
        mu_bins = measurement[1]


        # select correct lines to work with
        if file == "eu":
            #only use certain elements
            mu_bins = np.append(mu_bins[4:7],mu_bins[8:])
        peak_lines = np.array([mu_bins,ref_vals])
        return peak_lines[0], peak_lines[1]

    if detector == "ge":
        try:
            measurement = std.util.load_csv("../figs/"+file+"_ge_bin.csv", skiprows=1)
        except FileNotFoundError:
            return print("invalid file name!")
        mu_bins = measurement[1]

        if file == "eu":
            ids = [0, 1, 3, 16, 18, 19, 21, 23, 26]
            #print("hello", mu_bins[ids])
            mu_bins = mu_bins[ids]

        literature = lit_values_ge(str(file))

        ref_vals, ref_ints = literature
        #print(type(ref_errs))
        ref_vals = np.array(ref_vals)
        ref_ints = np.array(ref_ints)

        peak_lines = np.array([mu_bins,ref_vals])

        return peak_lines[0], peak_lines[1]


def lit_values_nai(element):
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


def lit_values_ge(element):
    cobalt = [1173.2,1332.5]
    caesium = [661.7]

    #europium = [121.7817]
    eu = std.util.load_csv("../figs/europium_lit.csv", skiprows=1)
    eu_energies = eu[0]
    eu_ints = eu[1]

    if str(element) == "co":
        return (cobalt, np.linspace(0,10,len(cobalt)))
    if str(element) == "cs":
        return (caesium, np.linspace(0,10,len(caesium)))
    if str(element) == "eu":
        return (eu_energies, eu_ints)
    return "error with the lit vals!"


def plot_data_ge(save_to=False):
    #get data for all elements
    x_meas_eu, y_lit_eu = get_data("eu", detector="ge")
    x_meas_co, y_lit_co = get_data("co", detector="ge")
    x_meas_cs, y_lit_cs = get_data("cs", detector="ge")

    #handle data and literature values
    x_temp = np.append(x_meas_co, x_meas_cs)
    x_meas = np.append(x_temp, x_meas_eu)
    x_meas_vals, x_meas_errs = p.ve(x_meas)

    #print("y lit:", y_lit_eu)

    # set placeholder errors to 0 for cs and co
    y_lit_co_errs = np.zeros_like(y_lit_co)
    y_lit_cs_errs = np.zeros_like(y_lit_cs)

    for i in range(len(y_lit_co)):
        y_lit_co[i] = p.ev(y_lit_co[i], y_lit_co_errs[i])
    for i in range(len(y_lit_cs)):
        y_lit_cs[i] = p.ev(y_lit_cs[i], y_lit_cs_errs[i])


    y_temp = np.append(y_lit_co, y_lit_cs)
    y_lit = np.append(y_temp, y_lit_eu)
    #print("y lit values:", y_lit)
    _, y_lit_errs = p.ve(y_lit)
    #y_lit_errs = np.append(y_lit_co_errs, y_lit_cs_errs)
    #y_lit_errs = np.append(y_lit_errs, y_lit_eu_err)
    #print(y_lit_errs)

    #fit linear function to data
    fit, cov = curve_fit(lambda x, a, b: a*x+b, x_meas_vals, ~y_lit*1e3)#, bounds=([0,0],[1000,1e10])) doesnt even need bounds this time!
    err = np.sqrt(np.diag(cov))
    #print(fit)

    #printable output
    params = []
    for i in range(len(fit)):
        params.append(p.ev(fit[i],err[i]))

    goodness = std.goodness_of_fit(~y_lit*1e3, std.linear(x_meas_vals,*fit))
    #print(goodness)
    x_ax = np.linspace(0, max(~x_meas), 500)

    #plot data and fit
    plt.errorbar(x_meas_vals, ~y_lit * 1e3, xerr=x_meas_errs, yerr=y_lit_errs * 1e3, label="zugeordnete Linien (HPGe-Detektor)",color="tab:red", **std.default.error_bar_def)
    plt.plot(x_ax,std.linear(x_ax,*fit), label=f"Anpassungsgerade, $R^2$={round(goodness, 5)}", color="tab:blue")

    plt.legend(loc="best", fontsize="small")
    std.default.plt_pretty("Kanalnr. / [Anzahl]","Energie / eV")

    print(params[0].format(),params[1].format())

    if save_to:
        plt.savefig("../figs/ge_calibration.pdf")

    else:
        plt.show()



    return

def plot_data_nai(save_to=False):

    x_meas_co, y_lit_co = get_data("co", detector="nai")
    x_meas_cs, y_lit_cs = get_data("cs", detector="nai")
    x_meas_eu, y_lit_eu = get_data("eu", detector="nai")

    x_temp = np.append(x_meas_co,x_meas_cs[:-1])
    x_meas = np.append(x_temp, x_meas_eu)
    x_meas_vals, x_meas_errs = p.ve(x_meas)
    y_temp = np.append(y_lit_co, y_lit_cs[:-1])
    y_lit = np.append(y_temp, y_lit_eu)
    print(y_lit)
    y_lit_errs = [0, 0, 0, 3e-4, 8e-4, 12e-4, 24e-4, 5e-3, 10e-3, 3e-3]

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
    goodness = std.goodness_of_fit(y_lit*1e3, std.linear(x_meas_vals,*fit))
    #print(fit)
    #print(err)
    x_ax = np.linspace(0, max(~x_meas), 500)
    #print(fit[0]*90+fit[1])

    plt.errorbar(x_meas_vals, y_lit*1e3, xerr=x_meas_errs, yerr=y_lit_errs, label="zugeordnete Linien (NaI-Detektor)",color="tab:red", **std.default.error_bar_def)
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

    plot_data_nai(save_to=False)

    return

if __name__ == "__main__":
    main()
