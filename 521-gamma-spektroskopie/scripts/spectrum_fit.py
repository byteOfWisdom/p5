# from scipy import curve_fit
import numpy as np
import std
from matplotlib import pyplot as plt
import copy


# def get_area(amplitudes):
#     amplitudes = np.convolve(amplitudes, np.ones(10), mode="same")
#     # plt.cla()
#     # plt.plot(amplitudes)
#     # plt.show()
#     max_id = np.where(amplitudes == max(amplitudes))[0][0]
#     area_size = 1
#     violation_counter = 0
#     while 1:
#         avg_amp = np.average([amplitudes[max_id + area_size], amplitudes[max_id - area_size]])
#         prev_avg_amp = np.average([amplitudes[max_id + (area_size - 1)], amplitudes[max_id - (area_size - 1)]])
#         if avg_amp >= prev_avg_amp:
#             violation_counter += 1
#         if avg_amp < 0.5 * amplitudes[max_id]:
#             violation_counter += 1
#         if violation_counter > 10:
#             break
#         area_size += 1
#     return area_size


# def get_area(amplitudes, chunk = 30):
#     # amplitudes = np.convolve(amplitudes, np.ones(30), mode="same")
    # max_id = np.where(amplitudes == max(amplitudes))[0][0]
#     x = max_id + chunk
#     while np.average(amplitudes[x - chunk:x]) > np.average(amplitudes[x:x + chunk]):
#         x += 1

#     y = max_id - chunk
#     while np.average(amplitudes[y:y + chunk]) > np.average(amplitudes[y - chunk:y]):
#         y -= 1

#     return y, x


def get_area(amplitudes):
    amplitudes = np.convolve(amplitudes, np.ones(10), mode="same")
    max_id = np.where(amplitudes == max(amplitudes))[0][0]
    is_big = amplitudes > 0.5 * np.average(amplitudes[max_id - 10: max_id + 10])
    upper = max_id
    while is_big[upper]:
        upper += 1
    lower = max_id
    while is_big[lower]:
        lower -= 1
    return lower, upper 


def fit_biggest_peak(channel, amplitude):
    max_id = np.where(amplitude == max(amplitude))[0][0]
    print(max_id)
    lower, upper = get_area(amplitude)
    p0 = [amplitude[max_id], channel[max_id], channel[max_id + upper] - channel[max_id]]
    channel_cut = channel[lower:upper]
    amplitude_cut = amplitude[lower:upper]
    fit_res, (_, goodness) = std.fit_func(std.gaussian, channel_cut, amplitude_cut, p0=p0, force_cf=True)
    plt.plot(channel_cut, amplitude_cut)
    return fit_res, goodness


def decomp_spectrum(channel, amplitude):
    fits = []
    goodness = []
    reduced_amps = copy.copy(amplitude)
    while max(reduced_amps) > 0.1 * max(amplitude):
        plt.plot(channel, reduced_amps)
        res, r_sq = fit_biggest_peak(channel, reduced_amps)
        print(r_sq)
        plt.plot(channel, std.gaussian(channel, *res), linestyle="dashed")
        plt.show()

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
        if goodness < 0.2:
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
    print(valid_lines)
    all_lines = std.make_n_gaussian(lc)
    res, _ = std.fit_func(all_lines, channel, amplitude, p0=valid_lines, force_cf=True)

    plt.plot(channel, amplitude)
    plt.plot(channel, all_lines(channel, *res))
    plt.show()
    return res
