from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import std
from sys import argv
import propeller as p

def get_data(file):
    data = np.transpose(np.loadtxt("../data/scope_"+file+".csv",skiprows=2,delimiter=","))
    return data[0], data[1], file

def rel_volts(volts):
    min_volt = min(volts)
    #print(min_volt)
    # corr_volts = []
    # for i in range(len(volts)):
    #     corr_volts.append(volts[i]+abs(min_volt))
    corr_volts = volts-min(volts)
    #print(corr_volts[50])
    return corr_volts

def plot_as_seperate(times, voltages,file,index,corrected=False):
    #plt.style.use("dark_background")
    # optionally correct voltages for zero position
    if corrected:
        voltages = rel_volts(voltages)
    plt.plot(times*1e6,voltages,color="tab:green",label="Aufnahme "+index,linewidth=0.8)
    std.default.plt_pretty(r"Zeit / $\mu$s","Spannung / V")
    plt.legend()
    #plt.savefig("../figs/oszi_uncorr_"+file+".pdf")
    plt.show()
    plt.close()
    return

def plot_as_subfigs(files,corrected=False):
    fig, axs = plt.subplots(nrows=2, ncols=2)
    for i in range(0,len(files)):
        #assign position on 2x2 grid
        # if i == 2:
        #     show = [5]
        first_ind = 0 if i in [1,2] else 1
        second_ind = 0 if i in [1,3] else 1

        times, volts, index = get_data(files[i])
        #optionally correct voltages for zero position
        if corrected:
            volts=rel_volts(volts)

        axs[first_ind,second_ind].plot(times*1e6,volts, color="tab:green",linewidth=0.8)
        print("plotting "+files[i])
        rec = str(i) if i != 0 else "4"
        axs[first_ind,second_ind].set_title("Aufnahme "+rec)
        axs[first_ind,second_ind].grid(which="major")
        axs[first_ind,second_ind].grid(which="minor", linestyle=":", linewidth=0.5)
        axs[first_ind,second_ind].minorticks_on()
        axs[first_ind,second_ind].set_xlabel(r"Zeit / $\mu$s")
        axs[first_ind,second_ind].set_ylabel("Spannung / V")
    fig.tight_layout()
    return

def plotting_figs(files,seperate=True):
    if seperate:
        runn_ind = 1
        for i in files:
            #print(runn_ind)
            index = plot_as_seperate(*get_data(i),corrected=False,index=str(runn_ind))
            runn_ind += 1
        return
    else:
        plot_as_subfigs(files)
        plt.savefig("../figs/oszi_uncorr_files_"+files[0]+"-"+files[-1]+".pdf")
        return

def analyse_shape(file):
    times, volts, _ = get_data(file)
    #print("file ",file,":")
    #plot_as_seperate(times, volts, argv[1], argv[1])

    times_div = 0.02e-6
    times_w_err = p.ev(times, times_div)
    volts_w_err = p.ev(volts, volts*0.05) # voltage err: 5% of value
    volts_corr = volts_w_err-min(volts_w_err)
    max_volt = max(volts_corr)
    max_time_index = np.where(volts_corr == max_volt)
    max_time = times_w_err[max_time_index]
    #print("maximum at:", max_time[0].format(), " index: ", max_time_index[0][0],"\n with: ", max_volt.format(),"V")

    sliced_volts = volts_corr[max_time_index[0][0]:]
    half_volt = max_volt/np.exp(1)
    half_rel_index = np.where(~sliced_volts < ~half_volt)[0][0]
    half_index = max_time_index+half_rel_index
    half_volt = volts_corr[half_index]
    half_time = times_w_err[half_index]
    #print("half height of ", half_volt[0][0].format(), "V at ", half_time[0][0].format(), "s" )

    duration = half_time - max_time
    #print("duration: ", duration[0][0].format(), "s")
    return max_volt, duration[0][0], max_time_index[0][0],half_index[0][0]

def av_shape(files):
    max_volts = []
    durations = []
    for i in files:
        max_volt, duration,_,_ = analyse_shape(i)
        max_volts.append(max_volt)
        durations.append(duration)
    av_max = sum(max_volts)/len(max_volts)
    av_duration = sum(durations)/len(durations)
    print("average height: ",av_max.format(), "V\naverage duration: ",av_duration.format(),"s")
    return

def exponential(x,a,b,c, t0):
    return c*np.exp(a*(x - t0))+b


def fit_shape(file):
    times, volts, _ = get_data(file)
    max_volt,duration, max_time_index, half_index = analyse_shape(file)
    init_guess = [-1e4,0.6,0.08,times[max_time_index]]
    upper_lim = int(half_index + (half_index - max_time_index))
    lower_lim = int(max_time_index-(max_time_index/8))

    #fit data to exponential
    fit, (sd, rsq) = std.fit_func(exponential, times[max_time_index:upper_lim], volts[max_time_index:upper_lim], p0=init_guess)
    params = p.ev(fit,sd)

    #data viz
    print("fit (c*np.exp(a*(x - d))+b):")
    for i in range(len(params)):
        print(params[i].format())
    print("R^2:",round(rsq,3))
    plt.plot(times[lower_lim:upper_lim+1300], volts[lower_lim:upper_lim+1300], label="gekürzte Messdaten (HPGe-Detektor)",color="orange")
    plt.plot(times[max_time_index:upper_lim], exponential(times[max_time_index:upper_lim], *fit),label=f"Anpassungsfunktion, $R^2$={round(rsq,3)}",color="darkblue")
    plt.legend(loc="best",fontsize="small")
    std.default.plt_pretty(r"Zeit / s","Spannung / V")
    #plt.show()
    plt.savefig("../figs/peak_curve_HPGe_"+file+".pdf")
    return



def main():
    #plotting_figs(argv[1:5],seperate=False)
    #correct files:
        # 6 7 8 9
        # 1 3 4 5
    # for i in range(1,len(argv)):
    #     analyse_shape(argv[i])
    #av_shape(argv[1:])
    fit_shape(argv[1])

    return

if __name__ == "__main__":
    main()
