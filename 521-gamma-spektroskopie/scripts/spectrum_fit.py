import numpy as np
import std
from matplotlib import pyplot as plt
import copy
import propeller as p
import scipy


plot_subfits = False


def get_area(amplitudes):
    max_id = np.where(amplitudes == max(amplitudes))[0][0]
    return get_area_around(amplitudes, max_id)


def get_area_around(amplitudes, index):
    is_big = amplitudes > 0.5 * amplitudes[index + 1]
    upper = index
    noise_rejection = 2
    while np.any(is_big[upper:upper + noise_rejection]):
        upper += 1
        if upper + noise_rejection >= len(amplitudes):
            break
    lower = index
    while np.any(is_big[lower - noise_rejection:lower]):
        lower -= 1
        if upper - noise_rejection >= 0:
            break
    return lower, upper


def fit_biggest_peak(channel, amplitude):
    max_id = np.where(amplitude == max(amplitude))[0][0]
    lower, upper = get_area(amplitude)
    p0 = [(channel[upper] - channel[lower]) / 2.5]
    channel_cut = channel[lower:upper]
    amplitude_cut = amplitude[lower:upper]
    fit_res, (_, goodness) = std.fit_func(
        lambda x, sigma: std.gaussian(x, amplitude[max_id], channel[max_id], sigma),
        channel_cut, amplitude_cut, p0=p0, force_cf=True)
    if plot_subfits:
        plt.plot(channel_cut, amplitude_cut, color="red")
    snr = np.log(amplitude[max_id]) * (upper - lower)
    return np.array([amplitude[max_id], channel[max_id], fit_res[0]]), snr


def strip_spectrum(channel, amplitude, ref_level):
    fits = []
    snr = []
    count = 0

    const = len(amplitude) // 1000
    ref_level = np.convolve(amplitude, np.ones(const), "same")

    definite_peaks = std.diff_find_maxima(amplitude, const)
    for i in range(len(definite_peaks)):
        rough_peak = definite_peaks[i]
        around = amplitude[rough_peak - const:rough_peak + const]
        peak = np.where(around == max(around))[0][0] + rough_peak - const
        definite_peaks[i] = peak

    plt.scatter(definite_peaks, amplitude[definite_peaks], color="purple")

    print(definite_peaks)
    for peak in definite_peaks:
        count += 1
        lower, upper = get_area_around(amplitude, peak)
        if upper - lower < 5:
            continue
        print(lower, upper)
        channel_cut = channel[lower:upper]
        amplitude_cut = amplitude[lower:upper]
        p0 = [(channel[upper] - channel[lower]) / 2.5]
        res, (_, r_sq) = std.fit_func(
            lambda x, sigma: std.gaussian(x, amplitude[peak], channel[peak], sigma),
            channel_cut, amplitude_cut, p0=p0, force_cf=True)
        res = [amplitude[peak], channel[peak], res[0]]
        if plot_subfits:
            plt.plot(channel, amplitude)
        if plot_subfits:
            print(r_sq)
            plt.plot(channel, std.gaussian(channel, *res), linestyle="dashed")
            plt.show()

        amplitude -= std.gaussian(channel, *res)
        amplitude[amplitude < 0] = 0
        fits.append(np.abs(res))
        snr.append(r_sq)

    while max(amplitude) > ref_level and count < 10:
        count += 1
        if plot_subfits:
            plt.plot(channel, amplitude)
        res, r_sq = fit_biggest_peak(channel, amplitude)
        if plot_subfits:
            print(r_sq)
            plt.plot(channel, std.gaussian(channel, *res), linestyle="dashed")
            plt.show()

        amplitude -= std.gaussian(channel, *res)
        amplitude[amplitude < 0] = 0
        fits.append(np.abs(res))
        snr.append(r_sq)
    return fits, snr


def make_spectrum_function(linecount, underground_fn):
    # all_lines = std.make_n_area_gaussian(linecount)
    all_lines = std.make_n_gaussian(linecount)
    return lambda x, *args: all_lines(x, *args[:3 * linecount]) + underground_fn(x, *args[3 * linecount:])


def plot_results(channel, amplitude, total_spectrum, underground_fn, res, ug_arg_count,  save_fig, goodness, lc):
    std.default.plt_pretty("Kanal", "Anzahl")
    plt.plot(channel, amplitude)
    plt.plot(channel, total_spectrum(channel, *res), label=f"$R^2={round(goodness, 3)}$")
    plt.legend()
    for i in range(lc):
        plt.plot(
            channel,
            std.gaussian(channel, res[3 * i], res[3 * i + 1], res[3 * i + 2]) + underground_fn(channel, *res[-ug_arg_count:]),
            linestyle="dashed")
    plt.vlines(res[1:-ug_arg_count:3], 0, max(amplitude), colors="green", linewidth=1)
    if save_fig:
        plt.savefig(save_fig)
        plt.cla()
    else:
        plt.show()


def decomp_spectrum(channel, amplitude, underground_fn, ug_arg_count, save_fig=False):
    fits = []
    snr = []

    reduced_snr = (np.average(amplitude) / np.sqrt(np.var(amplitude))) / 5
    print(reduced_snr)
    ug_fit, _ = std.fit_func(
        underground_fn, channel[amplitude < 0.75 * reduced_snr * max(amplitude)],
        amplitude[amplitude < 0.75 * reduced_snr * max(amplitude)], force_cf=True)


    reduced_amps = copy.copy(amplitude)
    reduced_amps = reduced_amps - underground_fn(channel, *ug_fit)
    fits, snr = strip_spectrum(channel, reduced_amps, reduced_snr * max(amplitude))

    # filter sensible peaks
    valid_lines = []
    lc = 0

    for params, snr in zip(fits, snr):
        # params[0] *= np.sqrt(2 * np.pi) * params[2]
        if snr < 0.1:
            print("this should not be printed!!")
            # continue
        if params[0] > 1.2 * max(amplitude):
            continue
        if params[1] > max(channel) or params[1] < min(channel):
            continue
        if params[2] < 1:
            continue
        valid_lines = np.append(valid_lines, params)
        lc += 1

    if plot_subfits:
        print(lc)
    total_spectrum = make_spectrum_function(lc, underground_fn)
    res, (std_devs, goodness) = std.fit_func(
        total_spectrum, channel, amplitude, y_errors=np.sqrt(amplitude),
        p0=np.append(valid_lines, np.zeros(ug_arg_count)), force_cf=True)

    plot_results(channel, amplitude, total_spectrum, underground_fn, res, ug_arg_count, save_fig, goodness, lc)

    res = p.ev(res, std_devs)
    res_dict = {
        "amp": res[0:-ug_arg_count:3],
        "mu": res[1:-ug_arg_count:3],
        "sigma": res[2:-ug_arg_count:3],
        "ug_params": res[-ug_arg_count:],
    }

    return res_dict, goodness
