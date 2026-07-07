#!python3
import numpy as np
import std
import sys

def path(n):
    if n < 10:
        n = f"0{n}"
    return f"../data/oszi/TRC{n}.CSV"


ch1_color = "tab:orange"
ch2_color = "tab:blue"
trace_alpha = 0.4
ax1 = std.default.plt.gca()
ax2 = ax1.twinx()

unit, factor = None, None


def plot_trace(n):
    global unit, factor
    trace_data = std.load_csv(path(n), delimiter=",", skiprows=1)
    t, ch1, ch2 = tuple(trace_data)
    if std.none(factor):
        unit, factor = find_time_unit(t)
    ax1.plot(t * factor, ch1, color=ch1_color, alpha=trace_alpha)
    ax2.plot(t * factor, ch2, color=ch2_color, alpha=trace_alpha)
    # ax1.scatter(t, ch1, color=ch1_color, alpha=trace_alpha, marker=".")
    # ax2.scatter(t, ch2, color=ch2_color, alpha=trace_alpha, marker=".")
    #print("finished", n)


def find_time_unit(times):
    t_scale = times[-1] - times[0]
    if t_scale <= 1000 * 1e-9:
        return "$ns$", 1e9
    if t_scale <= 1000 * 1e-6:
        return "$\\mu s$", 1e6
    if t_scale <= 1000 * 1e-3:
        return "$ms$", 1e3
    return "$s$", 1


def plot():
    global unit
    global trace_alpha
    ids = []
    out_file = None
    args = list(sys.argv)
    if "." in args[-1]:
        out_file = args[-1]
        args = args[:-1]

    for a in args[1:]:
        if "-" in a:
            lower = int(a.split("-")[0])
            upper = int(a.split("-")[1])
            ids += list(range(lower, upper+1))
            print("plotting ids:", ids)
        else:
            ids += [int(a)]
            trace_alpha = 1.0

    for id in ids:
        plot_trace(id)

    #print("done plotting all ids")
    ax1.set_xlabel(f"Zeit / {unit}")
    ax1.set_ylabel("Kanal 1 / V", color=ch1_color)

    ax1.grid(which="major")
    ax1.grid(which="minor", linestyle=":", linewidth=0.5)
    ax1.minorticks_on()
    ax2.set_ylabel("Kanal 2 / V", color=ch2_color)
    ax2.tick_params(axis="y",labelcolor=ch2_color)
    ax1.tick_params(axis="y",labelcolor=ch1_color)
    std.default.plt.gcf().set_size_inches(16/1.75, 9/1.75)

    if out_file:
        std.default.plt.savefig(out_file)
    else:
        std.default.plt.show()


if __name__ == "__main__":
    plot()
