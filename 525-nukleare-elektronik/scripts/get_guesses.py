import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import scipy
import std
import spectrum_fit
import propeller as p


class handle_keeper:
    def __init__(self):
        self.first_line = None
        self.second_line = None
        self.area = None
        self.curve = None

    def unrender(self):
        if self.first_line:
            self.first_line.remove()
        if self.second_line:
            self.second_line.remove()
        if self.area:
            self.area.remove()
        if self.curve:
            self.curve.remove()
        plt.draw()
    


class click_handler:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.interval = (min(y), max(y))
        self.lop = []
        self.handles = []
        self.total_handle = handle_keeper()
        self.p0s = []
        self.in_area = False
        self.total_params = None
        self.save = False


    def start(self):
        plt.connect("button_press_event", self.onclick)
        plt.connect("key_press_event", self.onkey)
        self.render()


    def render(self):
        plt.plot(self.x, self.y, linewidth=0.5)
        # plt.errorbar(self.x, self.y, np.sqrt(self.y), **std.default.error_bar_def)
        std.default.plt_pretty("Kanal", "Anzahl")
        plt.show()


    def fit_single_peak(self, a, b):
        a, b = min(a, b), max(a, b)
        curve = np.where((b >= self.x) & (self.x >= a))
        x_part, y_part = self.x[curve], self.y[curve]
        sigma_guess = (x_part[-1] - x_part[0]) / 2.5
        offset_guess = min(y_part)
        area_guess = max(y_part) * np.sqrt(2 * np.pi) * sigma_guess
        p0 = [area_guess, np.average(x_part), sigma_guess, offset_guess]
        try:
            res, cov = scipy.optimize.curve_fit(std.area_gaussian_ug, x_part, y_part, p0)
            red_chi_sq = std.reduced_chi_2(y_part, std.area_gaussian_ug(x_part, *res), res)#, sigma=np.sqrt(y_part))
            err = np.sqrt(np.diag(cov))
            self.p0s[-1] = list(p.ev(np.abs(res), err)) + [red_chi_sq]
            self.handles[-1].curve = plt.plot(x_part, std.area_gaussian_ug(x_part, *res), color="lightgreen")[0]
        except Exception as e:
            print(e)
            self.delete_last()


    def fit_double_peak(self, a, b):
        a, b = min(a, b), max(a, b)
        curve = np.where((b >= self.x) & (self.x >= a))
        x_part, y_part = self.x[curve], self.y[curve]
        sigma_guess = (x_part[-1] - x_part[0]) / 5
        offset_guess = min(y_part)
        area_guess = max(y_part) * np.sqrt(2 * np.pi) * sigma_guess
        mu_step = (x_part[-1] - x_part[0]) / 3
        mu_one = x_part[0] + mu_step
        mu_two = x_part[0] + 2 * mu_step
        p0 = [area_guess, area_guess, mu_one, mu_two, sigma_guess, sigma_guess, offset_guess]
        try:
            res, cov = scipy.optimize.curve_fit(std.double_area_gaussian, x_part, y_part, p0)
            red_chi_sq = std.reduced_chi_2(y_part, std.double_area_gaussian(x_part, *res), res)#, sigma=np.sqrt(y_part))
            err = np.sqrt(np.diag(cov))
            self.p0s[-1] = list(p.ev(np.abs(res), err)) + [red_chi_sq]
            self.handles[-1].curve = plt.plot(x_part, std.double_area_gaussian(x_part, *res), color="lightgreen")[0]
        except Exception as e:
            print(e)
            self.delete_last()


    def fit_spectrum(self):
        func = spectrum_fit.make_spectrum_function(len(self.p0s), spectrum_fit.poly_4)
        total_p0 = []
        for line in self.p0s:
            total_p0 += [line[0], line[1], line[2]]
        lower_bound = np.array(total_p0) - 0.2 * np.array(total_p0)
        upper_bound = np.array(total_p0) + 0.2 * np.array(total_p0)
        lower_bound = np.append(lower_bound, np.array([-1e3] * 4))
        upper_bound = np.append(upper_bound, np.array([1e3] * 4))
        total_p0 += [0] * 4

        res, cov = scipy.optimize.curve_fit(
            # std.make_n_area_gaussian(len(lines)),
            func,
            self.x, self.y,
            p0=total_p0,
            bounds=(lower_bound, upper_bound),
            xtol=1e-2,
            ftol=1e-2
        )
        self.total_handle.curve = plt.plot(self.x, func(self.x, *res), color="yellow")[0]
        plt.draw()
        print("updated total")


    def onclick(self, event):
        pass

    def delete_last(self):
        _ = self.lop.pop(-1)
        _ = self.p0s.pop(-1)
        handles = self.handles.pop(-1)
        handles.unrender()
        self.in_area = False


    def onkey(self, event):
        if event.key == "d":
            self.delete_last()

        elif event.inaxes and event.key == "w" and not self.in_area:
            self.lop += [(event.xdata, 0)]
            self.p0s += [[0, 0, 0]]
            self.handles.append(handle_keeper())
            self.handles[-1].first_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="green", linewidth=0.5)
            self.in_area = True
            plt.draw()

        elif event.inaxes and event.key == "w" and self.in_area:
            self.lop[-1] = (self.lop[-1][0], event.xdata)
            self.in_area = False
            self.handles[-1].second_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="green", linewidth=0.5)
            self.handles[-1].area = plt.fill_between(np.linspace(*self.lop[-1]), max(self.y), alpha=0.1, color="green")
            self.fit_single_peak(*self.lop[-1])
            plt.draw()

        elif event.inaxes and event.key == "W" and self.in_area:
            self.lop[-1] = (self.lop[-1][0], event.xdata)
            self.in_area = False
            self.handles[-1].second_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="green", linewidth=0.5)
            self.handles[-1].area = plt.fill_between(np.linspace(*self.lop[-1]), max(self.y), alpha=0.1, color="green")
            self.fit_double_peak(*self.lop[-1])
            plt.draw()

        elif event.key == "e":
            print("rerunning fit?")
            self.fit_spectrum()

        elif event.key == "S":
            total_p0 = []
            for line in self.p0s:
                total_p0 += [line[0], line[1], line[2], line[3], line[4]]
            lines_data ={
                "Fläche": total_p0[0::5],
                "$\\mu$": total_p0[1::5],
                "$\\sigma$": total_p0[2::5],
                "C": total_p0[3::5],
                "$\\chi^2_\\text{red}$": total_p0[4::5]
            }
            std.print_tex_table(lines_data, self.save + ".table")
            std.print_csv_table(lines_data, self.save + ".csv")
            plt.savefig(self.save + ".pdf")
            print("done saving!!")


def let_user_click_peaks(x_values, y_values, output):
    handler = click_handler(x_values, y_values)
    if output:
        handler.save = output
    handler.start()
    # ok might have wanted to do that earlier... let's not talk about that


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    # underground = np.transpose(np.loadtxt(argv[2]))
    # data[1] = data[1] - underground[1]
    out = False
    if len(argv) > 3:
        out = argv[3]
    let_user_click_peaks(data[0], data[1], out)
    return None

if __name__ == "__main__":
    main()
