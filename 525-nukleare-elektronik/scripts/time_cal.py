import std
import numpy as np
import propeller as p


def load_data(sigma=False):
    fname = "../data/fits/prompt.csv"
    data = std.load_csv(fname, skiprows=1)
    if sigma:
        return data[2]
    return data[1]

def time_cal(plot=False):
    lines = load_data()
    times = 16 * np.arange(5) + 8
    # res, (err, rsq) = std.curve_fit(std.linear, times, lines)
    res, (err, rsq) = std.odr_fit(std.linear, lines,  times)
    if plot:
        # std.default.plt_errorbar(times, lines)
        std.default.plt_errorbar(lines, times,marker="x")
        std.default.plt_func(std.linear, res, f"$R^2 = {round(rsq, 4)}$")
        # std.default.plt_finish("Zeit / ns", "Linie / Bins")
        std.default.plt_finish("Linie / Bins", "Zeit / ns", "../figs/time_cal.pdf")
        # std.default.plt_finish("Linie / Bins", "Zeit / ns")
    return p.ev(res, err)


def invert(a, b):
    return a, b
    # return 1 / a, -b / a


if __name__ == "__main__":
    res = time_cal(True)
    print(res)

    # std.write_file("../data/fits/time_cal_fit.txt", "$a = \\SI{" + res[0].format() + "}{\\nano\\second\\per\\bin}$ und $b = \\SI{" + res[1].format() + "}{\\nano\\second}$")
    std.write_file("../data/fits/time_cal_fit.txt", std.si_string("a", res[0], "\\nano\\second\\per\\bin") + " und " + std.si_string("b", res[1], "\\nano\\second"))

    sigmas = load_data(True)
    fwhm = np.sqrt(8 * np.log(2)) * sigmas
    time_res = fwhm * res[0]
    print(fwhm)
    print(time_res)
    table = {
        "\\sigma / Bin": sigmas,
        "FWHM / Bin": fwhm,
        "FWHM / ns": time_res
    }

    std.print_tex_table(table, "../latex/time_resolution.table")

    std.write_file("../data/fits/fwhm_time_res.txt", std.si_string("\\text{FWHM}" ,(sum(time_res) / len(time_res)).format(), "\\nano\\second"))
    print(sum(time_res) / len(time_res))
    
