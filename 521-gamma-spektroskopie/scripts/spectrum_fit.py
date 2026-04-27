import numpy as np
import std
from matplotlib import pyplot as plt
import copy


plot_subfits = False


def get_area(amplitudes):
    max_id = np.where(amplitudes == max(amplitudes))[0][0]
    is_big = amplitudes > 0.5 * amplitudes[max_id + 1]
    upper = max_id
    noise_rejection = 2
    while np.any(is_big[upper:upper + noise_rejection]):
        upper += 1
    lower = max_id
    while np.any(is_big[lower - noise_rejection:lower]):
        lower -= 1
    return lower, upper


def fit_biggest_peak(channel, amplitude):
    max_id = np.where(amplitude == max(amplitude))[0][0]
    lower, upper = get_area(amplitude)
    p0 = [(channel[upper] - channel[lower]) / 2.5]
    channel_cut = channel[lower:upper]
    amplitude_cut = amplitude[lower:upper]
    fit_res, (_, goodness) = std.fit_func(lambda x, sigma: std.gaussian(x, amplitude[max_id], channel[max_id], sigma), channel_cut, amplitude_cut, p0=p0, force_cf=True)
    if plot_subfits:
        plt.plot(channel_cut, amplitude_cut, color="red")
    snr = np.log(amplitude[max_id]) * (upper - lower)
    return np.array([amplitude[max_id], channel[max_id], fit_res[0]]), snr


def decomp_spectrum(channel, amplitude, underground_fn):
    fits = []
    snr = []

    ug_fit, _ = std.fit_func(underground_fn, channel[amplitude < 0.05 * max(amplitude)], amplitude[amplitude < 0.05 * max(amplitude)], force_cf=True)
    
    reduced_amps = copy.copy(amplitude)
    reduced_amps = reduced_amps - underground_fn(channel, *ug_fit)
    while max(reduced_amps) > 0.1 * max(amplitude):
        if plot_subfits:
            plt.plot(channel, reduced_amps)
        res, r_sq = fit_biggest_peak(channel, reduced_amps)
        if plot_subfits:
            print(r_sq)
            plt.plot(channel, std.gaussian(channel, *res), linestyle="dashed")
            plt.show()

        reduced_amps -= std.gaussian(channel, *res)
        reduced_amps[reduced_amps < 0] = 0
        fits.append(np.abs(res))
        snr.append(r_sq)


    print(max(amplitude))
    print(fits)
    # filter sensible peaks
    valid_lines = []
    lc = 0

    print("a    \t   mu  \t     sigma")
    for params, snr in zip(fits, snr):
        print(f"{params[0]}\t{params[1]}\t{params[2]}\t")
        if snr < 0.1:
            continue
        if params[0] > 1.2 * max(amplitude):
            continue
        if params[1] > max(channel) or params[1] < min(channel):
            continue
        if params[2] < 1:
            continue
        valid_lines = np.append(valid_lines, params)
        lc += 1

    print(lc)
    all_lines = std.make_n_gaussian(lc)
    res, _ = std.fit_func(all_lines, channel, amplitude, y_errors=np.sqrt(amplitude), p0=valid_lines, force_cf=True)

    plt.plot(channel, amplitude)
    plt.plot(channel, all_lines(channel, *res))
    plt.vlines(res[1::3], 0, 500, colors="yellow", linestyles="dashed")
    plt.show()
    return res
