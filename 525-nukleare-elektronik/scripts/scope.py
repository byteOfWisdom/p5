#!/usr/bin/python3
import numpy as np
from sys import argv
import std

def path(n):
    if n < 10:
        n = f"0{n}"
    return f"../data/oszi/TRC{n}.CSV"


ch1_color = "tab:orange"
ch2_color = "tab:blue"
trace_alpha = 0.4
ax1 = std.default.plt.gca()
ax2 = ax1.twinx()

def plot_trace(n):
    trace_data = std.load_csv(path(n), delimiter=",", skiprows=1)
    t, ch1, ch2 = tuple(trace_data)
    ax1.plot(t, ch1, color=ch1_color, alpha=trace_alpha)
    ax2.plot(t, ch2, color=ch2_color, alpha=trace_alpha)
    # ax1.scatter(t, ch1, color=ch1_color, alpha=trace_alpha, marker=".")
    # ax2.scatter(t, ch2, color=ch2_color, alpha=trace_alpha, marker=".")


if __name__ == "__main__":
    ids = []
    for a in argv[1:]:
        if "-" in a:
            lower = int(a.split("-")[0])
            upper = int(a.split("-")[1])
            ids += list(range(lower, upper))
        else:
            ids += [int(a)]

    for id in ids:
        plot_trace(id)
    ax1.set_xlabel("Zeit / s")
    ax1.set_ylabel("Kanal 1 / V", color=ch1_color)

    ax1.grid(which="major")
    ax1.grid(which="minor", linestyle=":", linewidth=0.5)
    ax1.minorticks_on()
    ax2.set_ylabel("Kanal 2 / V", color=ch2_color)
    ax2.tick_params(axis="y",labelcolor=ch2_color)
    ax1.tick_params(axis="y",labelcolor=ch1_color)
    std.default.plt.show()
    
