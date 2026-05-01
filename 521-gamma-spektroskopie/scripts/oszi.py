from matplotlib import pyplot as plt
import numpy as np
import std
from sys import argv

def get_data(file):
    data = np.transpose(np.loadtxt("../data/scope_"+file+".csv",skiprows=2,delimiter=","))
    return data[0], data[1], file

def rel_volts(volts):
    min_volt = min(volts)
    #print(min_volt)
    corr_volts = []
    for i in range(len(volts)):
        corr_volts.append(volts[i]+abs(min_volt))
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
    plt.savefig("../figs/oszi_uncorr_"+file+".pdf")
    plt.close()
    return


def plot_as_subfigs(files,corrected=False):
    fig, axs = plt.subplots(nrows=2, ncols=2)
    for i in range(0,len(files)):
        #assign position on 2x2 grid
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

def analyse_shape():
    # TO DO
    return


def main():
    plotting_figs(argv[1:5],seperate=False)
    #correct files:
        # 6 7 8 9
        # 1 3 4 5
    return

if __name__ == "__main__":
    main()
