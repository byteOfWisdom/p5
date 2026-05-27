import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import std
import propeller as p
from sys import argv

def get_data(file):
    times, currents = np.genfromtxt(file, delimiter=",",skip_header=18, usecols=(3,4),autostrip=True, unpack=True)
    #print(times[:10])
    #print(currents[:10])
    #currents.strip()
    return times, currents

def plot_data(file):
    times, currents = get_data(file)
    min_time, half_time, duration = analyze(file)
    #print(char_times)
    #print(times[:10])
    #print(currents[:10])
    t_err = []
    c_err = []
    for t in range(len(times)):
        t_err.append(2*1e-9) #set times err to time div/2
        c_err.append(abs(currents[t]*3e-2)) #set current errs to 3%


    plt.vlines([~min_time, ~half_time], min(currents-0.01), max(currents+0.01), color="tab:green", label=f"Dauer: {duration.format()} s")
    plt.axvspan(~min_time, ~half_time, ymin = 0, ymax = 1, alpha=0.3, color='lightgreen', linewidth=0)

    #plt.errorbar(times, currents, xerr=t_err, yerr=c_err, label="Messdaten", **std.default.error_bar_def)
    plt.plot(times, currents, label="Messdaten", color = "tab:orange")
    #rsq = fit_shape(file)
    std.default.plt_pretty("Zeit / s", "Spannung / V")
    plt.legend(loc="best", fontsize="small")
    plt.show()
    #plt.savefig("../figs/"+argv[1].strip("../data/").strip(".CSV")+".pdf")

    return

def fit_shape(file):
    times, currents = get_data(file)

    min_time, half_time, duration = analyze(file, verbose=False)

    min_time_index = np.where(times == ~min_time)[0][0]
    half_time_index = np.where(times == ~half_time)[0][0]

    init_guess = [1e4,0.6,-0.08,~min_time]
    upper_lim = int(half_time_index + (half_time_index - min_time_index))
    lower_lim = int(min_time_index-(min_time_index/8))

    #fit data to exponential
    fit, (sd, rsq) = std.curve_fit(exponential, times[min_time_index:upper_lim], currents[min_time_index:upper_lim], p0=init_guess)
    params = p.ev(fit,sd)

    #data viz
    print("fit (c*np.exp(a*(x - d))+b):")
    names = ["a", "b", "c", "t0"]
    for i in range(len(params)):
        print(names[i],": ", params[i].format())
    print("R^2:",round(rsq,3))
    plt.plot(times[min_time_index:upper_lim], exponential(times[min_time_index:upper_lim], *fit),label=f"Anpassungsfunktion, $R^2$={round(rsq,3)}",color="darkblue")
    return rsq

def exponential(x,a,b,c, t0):
    return c*np.exp(a*(x - t0))+b

def analyze(file, verbose=True):
    times, currents = get_data(file)
    t_err = []
    c_err = []
    for t in range(len(times)):
        t_err.append(2*1e-9) #set times err to time div/2
        c_err.append(abs(currents[t]*3e-2)) #set current errs to 3%

    times_w_err = p.ev(times, t_err)
    currents_w_err = p.ev(currents, c_err)

    min_volt = min(currents_w_err)
    half_volt_height = min_volt/np.exp(1)

    min_time_index = np.where(currents == ~min_volt)
    #print(min_time_index[0][0])
    min_time = times_w_err[min_time_index]

    sliced_curr = currents_w_err[min_time_index[0][0]:]

    half_rel_index = np.where(~sliced_curr >= ~half_volt_height)[0][0]
    #print(half_rel_index)
    half_index = min_time_index + half_rel_index
    #print(half_index[0][0])
    half_volt = currents_w_err[half_index]
    half_time = times_w_err[half_index]

    duration = half_time[0][0] - min_time[0]

    if verbose is True:
        print("min volt height: ", ~min_volt, "half volt height: ", ~half_volt_height) #sanity check
        print("min time: ", ~min_time[0],"half height time: ", ~half_time[0][0])
        print("duration: ", duration.format())

    return min_time[0], half_time[0][0], duration

def drift_signal(file):
    times, currents = get_data(file)
    t_err = []
    c_err = []
    for t in range(len(times)):
        t_err.append(2*1e-9) #set times err to time div/2
        c_err.append(abs(currents[t]*3e-2)) #set current errs to 3%

    times_w_err = p.ev(times, t_err)
    currents_w_err = p.ev(currents, c_err)

    max_curr = max(currents_w_err)

    signal_level = p.ev(-0.55,0.55*3e-2) #estimated as average of measurements before signal w/ usual 3% error on current
    #print(signal_level.format())

    start_where = np.where(~signal_level <= currents)
    x_start = times[517] # for file 004!!!
    x_start = p.ev(x_start, 4*1e-9) #estimated error of double the div
    #x_start = p.ev(-7.150e-9, 4*1e-9) #estimated error of double the div

    list_where = np.where(~signal_level == currents)[0]
    list_len = len(list_where)
    #print(list_where)

    if list_len != 0:
       x_end = p.ev(list_where[:-1],4*1e-9)
    else:
        x_end = p.ev(1.23184e-8, 4*1e-9)
    #print(x_end.format())
    duration = x_end - x_start
    print("duration:", duration.format(), "s")



    plt.plot(times, currents, label = "Messdaten", color="tab:blue")
    plt.hlines([~signal_level], min(times), max(times), label = f"geschätztes Signallevel: {signal_level.format()} V", color="tab:orange")
    #plt.hlines([~max_curr], min(times), max(times), label = "lowest amplitude of current", color="tab:green")
    plt.vlines([~x_start, ~x_end], 0, min(currents), color = "tab:green", label = f"Signaldauer: {duration.format()} s")
    plt.axvspan(~x_start, ~x_end, ymin = 0, ymax = 1, alpha=0.3, color='lightgreen', linewidth=0)
    #plt.axvspan(1, 2, ymin = 0, ymax = 1, alpha=0.3, color='lightgreen', linewidth=0)

    #plt.vlines(start_where, 0, min(currents), color = "tab:green")
    plt.legend(fontsize = "small", loc = "lower right")
    std.default.plt_pretty("Zeit / s", "Spannung / V")
    #plt.show()
    plt.savefig("../figs/"+file.strip(".csv").strip("../data/")+".pdf")


    return


def main():
    #plot_data(argv[1])
    drift_signal(argv[1])
    return

if __name__ == "__main__":
    main()
