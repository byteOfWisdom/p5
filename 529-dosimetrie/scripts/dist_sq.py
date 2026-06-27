# %%
import std
import tomllib
from matplotlib import pyplot as plt
import numpy as np
import propeller as p



def convert(x):
    return x / 32.4


def quadratic(x, a, c):
    denom = (x ** 2)
    return (a / denom) + c
# %%
def get_data(number: str):
    fhandle = open("../data/meta.toml", "rb")
    meta = tomllib.load(fhandle)
    fhandle.close()
    fname = "../data/" + meta["abstand_sq"]["messung_1"]["csv_file"] if number == "1" else "../data/" + meta["abstand_sq"]["messung_2"]["csv_file"]

    data = std.util.load_csv(fname, skiprows=1)

    distance = p.ev(data[3] + 11, 0.5)
    dose_before = p.ev(data[1], 0.1)
    dose_after = p.ev(data[2], 0.1)
    time = p.ev(data[0], 0.5)
    #print(f"i have {len(distance)} points")
    #print(f"done getting data for {number}")
    return distance, dose_before, dose_after, time, number

def handle_data(distance, dose_before, dose_after, time, number):
    dose_rate = (dose_after - dose_before) / time

    equiv_rate = convert(dose_rate)

    print(r"\hline")
    print(r"Abstand / cm & Dauer / s & $H_\text{1}$ / mSv & $H_\text{2}$ / mSv & h / mSv/s \\")
    print(r"\hline")
    for i in range(len(dose_rate)):
        print(r"\num{", distance[i].format(), r"} & \num{", time[i].format(), r"} & \num{", dose_before[i].format(), r"} & \num{", dose_after[i].format(), r"} & \num{", dose_rate[i].format(), r"} \\")
    print(r"\hline")

    fit, (errs, rsq) = std.curve_fit(quadratic, distance, dose_rate)
    params = p.ev(fit, errs)


    #print("parameter a & c:")
    #[print(params[i].format()) for i in range(len(fit))]

    messung = "Mo" if number == "1" else "Cu"

   # print(f"i have {len(dose_rate)} points for {messung}")
    plt.errorbar(~distance, ~dose_rate, xerr = p.ve(distance)[1], yerr = p.ve(dose_rate)[1], label = f"Messwerte {messung}",  **std.default.error_bar_def)
    plt.scatter(~distance, ~dose_rate)
    xrange = np.linspace(min(~distance) - 2, max(~distance) + 2)
    plt.plot(xrange, quadratic(xrange, *fit), linewidth = 0.9, label = rf"Anpassungsfunktion {messung}, $R^2$={round(rsq,3)}")
    return

def main():
    data2 = get_data("2")
    handle_data(*data2)
    #data1 = get_data("1")
    #handle_data(*data1)
    plt.legend()
    std.default.plt_pretty("Abstand / cm", "Dosisleistung / mSv/s")
    #plt.savefig("../figs/dist_sq.pdf")
    #plt.show()

if __name__ == "__main__":
    main()
