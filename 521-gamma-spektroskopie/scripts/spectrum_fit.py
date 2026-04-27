import numpy as np
import std
from matplotlib import pyplot as plt
import copy
import propeller as p

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


def decomp_spectrum(channel, amplitude, underground_fn, ug_arg_count, save_fig=False):
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


    # filter sensible peaks
    valid_lines = []
    lc = 0

    # print("a    \t   mu  \t     sigma")
    for params, snr in zip(fits, snr):
        # print(f"{params[0]}\t{params[1]}\t{params[2]}\t")
        if snr < 0.1:
            print("this should not be printed!!")
            continue
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
    all_lines = std.make_n_gaussian(lc)
    with_ug = lambda x, *args: all_lines(x, *args[:3 * lc]) + underground_fn(x, *args[3 * lc:])
    res, (std_devs, goodness) = std.fit_func(with_ug, channel, amplitude, y_errors=np.sqrt(amplitude), p0=np.append(valid_lines, np.zeros(ug_arg_count)), force_cf=True)

    std.default.plt_pretty("Kanal", "Anzahl")
    plt.plot(channel, amplitude)
    plt.plot(channel, with_ug(channel, *res), label=f"$R^2={round(goodness, 3)}$")
    plt.legend()
    for i in range(lc):
        plt.plot(channel, std.gaussian(channel, res[3 * i], res[3 * i + 1], res[3 * i + 2]) + underground_fn(channel, *res[-ug_arg_count:]), linestyle="dashed")
    plt.vlines(res[1:-ug_arg_count:3], 0, 500, colors="green", linestyles="dotted", linewidth=1)
    if save_fig:
        plt.savefig(save_fig)
        plt.cla()
    else:
        plt.show()

    res = p.ev(res, std_devs)
    res_dict = {
        "amp": res[0:-ug_arg_count:3],
        "mu": res[1:-ug_arg_count:3],
        "sigma": res[2:-ug_arg_count:3],
        "ug_params": res[-ug_arg_count:],
    }
    
    return res_dict, goodness
