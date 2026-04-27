# from scipy import curve_fit
import numpy as np
import std
from matplotlib import pyplot as plt
import copy


def get_area(amplitudes):
    amplitudes = np.convolve(amplitudes, np.ones(10), mode="same")
    max_id = np.where(amplitudes == max(amplitudes))[0][0]
    is_big = amplitudes > 0.5 * np.average(amplitudes[max_id - 10: max_id + 10])
    upper = max_id
    while is_big[upper] or is_big[upper + 1]:
        upper += 1
    lower = max_id
    while is_big[lower] or is_big[lower - 1]:
        lower -= 1
    return lower, upper 


def fit_biggest_peak(channel, amplitude):
    max_id = np.where(amplitude == max(amplitude))[0][0]
    lower, upper = get_area(amplitude)
    p0 = [amplitude[max_id], channel[max_id], channel[max_id + upper] - channel[max_id]]
    channel_cut = channel[lower:upper]
    amplitude_cut = amplitude[lower:upper]
    fit_res, (_, goodness) = std.fit_func(std.gaussian, channel_cut, amplitude_cut, p0=p0, force_cf=True)
    # plt.plot(channel_cut, amplitude_cut)
    return fit_res, goodness


def decomp_spectrum(channel, amplitude, underground_fn):
    fits = []
    goodness = []

    ug_fit, _ = std.fit_func(underground_fn, channel[amplitude < 0.05 * max(amplitude)], amplitude[amplitude < 0.05 * max(amplitude)], force_cf=True)
    
    reduced_amps = copy.copy(amplitude)
    reduced_amps = reduced_amps - underground_fn(channel, *ug_fit)
    plt.plot(channel, underground_fn(channel, *ug_fit))
    while max(reduced_amps) > 0.1 * max(amplitude):
        # plt.plot(channel, reduced_amps)
        res, r_sq = fit_biggest_peak(channel, reduced_amps)
        print(r_sq)
        # plt.plot(channel, std.gaussian(channel, *res), linestyle="dashed")
        # plt.show()

        reduced_amps -= std.gaussian(channel, *res)
        reduced_amps[reduced_amps < 0] = 0
        fits.append(np.abs(res))
        goodness.append(r_sq)


    print(max(amplitude))
    print(fits)
    # filter sensible peaks
    valid_lines = []
    lc = 0

    for params, goodness in zip(fits, goodness):
        # if goodness < 0.2:
        #     continue
        if params[0] > 1.2 * max(amplitude):
            continue
        if params[1] > max(channel) or params[1] < min(channel):
            continue
        if params[2] < 1:
            continue
        valid_lines = np.append(valid_lines, params)
        lc += 1

    print(lc)
    print(valid_lines)
    all_lines = std.make_n_gaussian(lc)
    res, _ = std.fit_func(all_lines, channel, amplitude, p0=valid_lines, force_cf=True)

    plt.plot(channel, amplitude)
    plt.plot(channel, all_lines(channel, *res))
    plt.vlines(res[1::3], 0, 500, colors="yellow", linestyles="dashed")
    plt.show()
    return res
