import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import scipy
import std
import spectrum_fit


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


    def start(self):
        plt.connect("button_press_event", self.onclick)
        plt.connect("key_press_event", self.onkey)
        self.render()


    def render(self):
        plt.plot(self.x, self.y)
        std.default.plt_pretty("Kanal", "Anzahl")
        plt.show()


    def fit_single_peak(self, a, b):
        a, b = min(a, b), max(a, b)
        curve = np.where((b >= self.x) & (self.x >= a))
        x_part, y_part = self.x[curve], self.y[curve]
        p0 = [max(y_part), np.average(x_part), (x_part[-1] - x_part[0]) / 2.5]
        try:
            res, _ = scipy.optimize.curve_fit(std.gaussian, x_part, y_part, p0)
            self.p0s[-1] = np.abs(res)
            self.handles[-1].curve = plt.plot(x_part, std.gaussian(x_part, *res), color="darkgreen")[0]
        except Exception as e:
            print(e)
            self.delete_last()


    def fit_spectrum(self):
        func = spectrum_fit.make_spectrum_function(len(self.p0s), spectrum_fit.poly_4)
        total_p0 = []
        for p in self.p0s:
            total_p0 += [p[0] * (p[2] * np.sqrt(2 * np.pi)), p[1], p[2]]
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
        # p0 += [0, 0, 0, 0]
        # res, (err, goodness) = std.fit_func(spectrum_fit.make_spectrum_function(len(self.lop), spectrum_fit.poly_4), self.x, self.y, y_errors=np.sqrt(self.y), p0=p0, force_cf=True)
        self.total_handle.curve = plt.plot(self.x, func(self.x, *res), color="lightgreen")[0]
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
            self.handles[-1].first_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="red")
            self.in_area = True
            plt.draw()

        elif event.inaxes and event.key == "w" and self.in_area:
            self.lop[-1] = (self.lop[-1][0], event.xdata)
            self.in_area = False
            self.handles[-1].second_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="green")
            self.handles[-1].area = plt.fill_between(np.linspace(*self.lop[-1]), max(self.y), alpha=0.25, color="green")
            self.fit_single_peak(*self.lop[-1])
            plt.draw()

        if event.key == "e":
            print("rerunning fit?")
            self.fit_spectrum()


def let_user_click_peaks(x_values, y_values):
    handler = click_handler(x_values, y_values)
    handler.start()
    # ok might have wanted to do that earlier... let's not talk about that


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    # out = argv[2]
    let_user_click_peaks(data[0], data[1])
    return None

if __name__ == "__main__":
    main()
