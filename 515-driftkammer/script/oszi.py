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


#doesnt work for drift chamber measurements due to 2nd & consec. peaks overlaying each other -> needs other method for determining half life length
def analyze_second(file):
    times, currents = get_data(file)
    t_err = []
    c_err = []
    for t in range(len(times)):
        t_err.append(2*1e-9) #set times err to time div/2
        c_err.append(abs(currents[t]*3e-2)) #set current errs to 3%
    times_w_err = p.ev(times, t_err)
    currents_w_err = p.ev(currents, c_err)

    min_time_1, half_time_1, duration_1 = analyze(file, verbose=False)
    print("min time 1: ", min_time_1.format())
    min_time_index_1 = np.where(times == ~min_time_1)[0][0]
    half_time_index_1 = np.where(times == ~half_time_1)[0][0]

    plt.plot(times, currents, label="Messdaten", color = "tab:orange")

    plt.vlines([~min_time_1, ~half_time_1], min(currents-0.01), max(currents+0.01), color="tab:green", label=f"Dauer: {duration_1.format()} s")
    plt.axvspan(~min_time_1, ~half_time_1, ymin = 0, ymax = 1, alpha=0.3, color='lightgreen', linewidth=0)


    sliced_times = times_w_err[half_time_index_1:]
    sliced_currents = currents_w_err[half_time_index_1:]

    min_volt_2 = min(sliced_currents)
    half_volt_height_2 = min_volt_2/np.exp(1)

    min_time_index_2 = np.where(currents_w_err == min_volt_2)
    print("min time index 2: ", min_time_index_2[0][0])
    min_time_2 = times_w_err[min_time_index_2[0][0]]
    print("min time 2: ", min_time_2.format())

    # sliced_curr_2 = currents_w_err[min_time_index_2[0][0]:]
    # half_rel_index_2 = np.where(~sliced_curr_2 >= ~half_volt_height_2)[0][0]
    # half_index_2 = min_time_index_2 + half_rel_index_2
    # half_volt_2 = currents_w_err[half_index_2]
    # half_time_2 = times_w_err[half_index_2]

    # duration = half_time_2[0][0] - min_time_2[0]



    # plt.vlines([~min_time_2, ~half_time_2], min(currents-0.01), max(currents+0.01), color="tab:green", label=f"Dauer: {duration_1.format()} s")
    # plt.axvspan(~min_time_2, ~half_time_2, ymin = 0, ymax = 1, alpha=0.3, color='lightgreen', linewidth=0)

    # std.default.plt_pretty("Zeit / s", "Spannung / V")
    plt.legend(loc="best", fontsize="small")
    plt.show()

    return



def main():
    #analyze_second(argv[1])
    #analyze(argv[1])
    plot_data(argv[1])
    return

if __name__ == "__main__":
    main()
