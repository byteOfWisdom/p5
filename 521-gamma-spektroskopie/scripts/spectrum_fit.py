import numpy as np
import std
from matplotlib import pyplot as plt
import propeller as p
import scipy

plot_subfits = False


def make_spectrum_function(linecount, underground_fn):
    # all_lines = std.make_n_area_gaussian(linecount)
    all_lines = std.make_n_area_gaussian(linecount)
    return lambda x, *args: all_lines(x, *args[:3 * linecount]) + underground_fn(x, *args[3 * linecount:])


def make_peak_view(amplitude):
    const = len(amplitude) // 100
    underground = np.convolve(amplitude, np.ones(15 * const) / (15 * const), "same")
    baseline_level = amplitude - underground
    res = amplitude - amplitude
    for sigma in np.linspace(0.1, 0.75, 20):
        kernel = std.gaussian(np.linspace(0, 2, len(amplitude)), 1, 1, sigma / const)
        baseline_level = np.convolve(baseline_level, kernel, "same") / const
        res +=baseline_level
    kernel = std.gaussian(np.linspace(0, 2, len(amplitude)), 1, 1, 0.05 / const)
    for _ in range(3):
        res = np.convolve(res, kernel, "same") / sum(kernel)
    return res * (1 / max(res))
    return baseline_level
    


def is_gaussian_peak(values):
    if len(values) < 3:
        return False, [0, 0, 0]
    p0 = [max(values), len(values) / 2, len(values) / 2.5]
    gauss_res, (_, gaussian_rsq) = std.fit_func(
                                                std.gaussian,
                                                np.arange(len(values)), values,
                                                y_errors=np.sqrt(values),
                                                force_cf=True, p0=p0)
    lin_res, (_, linear_rsq) = std.fit_func(
                                            std.linear,
                                            np.arange(len(values)), values,
                                            y_errors=np.sqrt(values),
                                            force_cf=True,
                                            p0=[values[0], (values[-1] - values[0]) / len(values)])
    # quad_res, (_, quad_rsq) = std.fit_func(std.quadratic, np.arange(len(values)), values, force_cf=True, p0=p0)
    # plt.plot(values)
    # plt.plot(std.gaussian(np.arange(len(values)), *gauss_res))
    # plt.plot(std.linear(np.arange(len(values)), *lin_res))
    # # plt.title(gaussian_rsq - linear_rsq)
    # plt.title(f"{gaussian_rsq}  {len(values) / gauss_res[2]}")
    # plt.show()
    is_gauss = len(values) / gauss_res[2] > 1. and gaussian_rsq > 0.3 and gaussian_rsq > linear_rsq
    return is_gauss, gauss_res


def find_gaussian_peaks(amplitude):
    peak_view = make_peak_view(amplitude)
    const = len(amplitude) // 1000
    definite_peaks, props = scipy.signal.find_peaks(peak_view, width=const // 2, prominence=1e-3)

    # print(props)
    gaussian_shaped, params = [], []

    for peak, width in zip(definite_peaks, props["widths"]):
        # print(peak, width)
        is_gauss, res = is_gaussian_peak(amplitude[peak - int(width * 1.5) // 2: peak + int(width * 1.5) // 2])
        if is_gauss:
            gaussian_shaped.append(peak)
            p0 = (res[0] * np.sqrt(np.pi * 2) * res[2], peak, res[2])
            params.append(p0)

    # plt.cla()
    # plt.plot(channel, amplitude, linewidth=0.2)
    # # plt.scatter(channel[definite_peaks], props["prominences"], color="red")
    # # plt.plot(channel, peak_view * max(amplitude))# / max(peak_view)))
    # plt.scatter(channel[gaussian_shaped], amplitude[gaussian_shaped], color="green")
    # plt.title("after first round gaussian fitting")
    # plt.show()

    return gaussian_shaped, params


def analyze_spectrum(x_values, y_values):
    lines, rough_params = find_gaussian_peaks(y_values)
    total_p0 = []
    for a, mu, sigma in rough_params:
        p0 = [a, x_values[mu], x_values[mu + int(sigma)] - x_values[mu]]
        total_p0 += p0

    lower_bound = 0.5 * np.array(total_p0 + [-np.inf] * 3)
    upper_bound = 2 * np.array(total_p0 + [np.inf] * 3)
    res, _ = scipy.optimize.curve_fit(
        std.make_n_area_gaussian(len(lines)),
        x_values, y_values,
        p0=total_p0,
        bounds=(lower_bound, upper_bound),
        xtol=1e-5,
        ftol=1e-5
    )

    plt.plot(x_values, y_values, linewidth=0.5)
    plt.plot(x_values, std.make_n_area_gaussian(len(lines))(x_values, *res), label="final")
    std.default.plt_finish("channel", "count")
    # return p.ev(res, err), rsq
