import numpy as np
import propeller as p
import std
from matplotlib import pyplot as plt

def get_upper_files(broken=False):
    upper_files = []
    upper_pos = [-6, 6, -12, 12, -24, 24, 0]
    position = 0
    pos_colors = ["tab:green", "tab:red"]

    if broken:
        print("for broken upper row")
        for i in [3,4]:
            data = np.transpose(np.loadtxt(f"../data/oszi/TRC0{str(i)}.CSV", delimiter=",", skiprows = 1))
            time = data[0]
            voltage = data[1]*10
            name = f"{upper_pos[position]}"
            plt.plot(time, voltage, label = f"{name} V")
            plt.hlines(int(upper_pos[position]), min(time)-1e-3, max(time)+1e-3, linestyle = "dashed", label = f"erwartet für {name} V", color = pos_colors[position])
            position += 1
    else:
        print("for working upper row")
        position = 2
        for i in range(5, 10):
            #print(i)
            data = np.transpose(np.loadtxt(f"../data/oszi/TRC0{str(i)}.CSV", delimiter=",", skiprows = 1))
            time = data[0]
            voltage = data[1]
            #print(position)
            name = "Erde" if upper_pos[position] == 0 else f"{upper_pos[position]} V"
            plt.plot(time, voltage, label = f"{name}")
            #plt.hlines(int(upper_pos[position]), min(time)-1e-3, max(time)+1e-3, linestyle = "dashed", color = "brown")
            position += 1

    #plotting figure
    plt.legend(fontsize = "small", loc = "best")
    std.default.plt_pretty("Zeit / s", "Spannung / V")
    if broken:
        plt.savefig("../figs/crates_broken.pdf")
    else:
        plt.savefig("../figs/crates_upper.pdf")

def get_lower_files():
    upper_files = []
    upper_pos = [24, -24, 12, -12, 6, -6, 0]
    position = 0
    print("for lower row")
    for i in range(10, 17):
        data = np.transpose(np.loadtxt(f"../data/oszi/TRC{str(i)}.CSV", delimiter=",", skiprows = 1))
        time = data[0]
        voltage = data[1]
        name = "Erde" if upper_pos[position] == 0 else f"{upper_pos[position]} V"
        plt.plot(time, voltage, label = f"{name}")
        position += 1
    plt.legend(fontsize = "small", loc = "best")
    std.default.plt_pretty("Zeit / s", "Spannung / V")
    plt.savefig("../figs/crates_lower.pdf")


def main():
    get_upper_files(broken=True)
    get_upper_files(broken=False)
    get_lower_files()
    print("done")
    return

if __name__ == "__main__":
    main()
