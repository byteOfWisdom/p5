import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import scipy
import std


class handle_keeper:
    def __init__(self):
        self.first_line = None
        self.second_line = None
        self.area = None

    def unrender(self):
        if self.first_line:
            self.first_line.remove()
        if self.second_line:
            self.second_line.remove()
        if self.area:
            self.area.remove()
        plt.draw()
    


class click_handler:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.interval = (min(y), max(y))
        self.lop = []
        self.handles = []
        self.in_area = False


    def start(self):
        plt.connect("button_press_event", self.onclick)
        plt.connect("key_press_event", self.onkey)
        self.render()


    def render(self):
        plt.plot(self.x, self.y)
        std.default.plt_pretty("Kanal", "Anzahl")
        plt.show()


    def fit_single_peak(self, a, b):
        curve = np.where(b >= self.x >= a)
        x_part, y_part = self.x[curve], self.y[curve]
        p0 = [max(y_part), np.average(x_part), (x_part[-1] - x_part[0]) / 2.5]
        # res, _ = scipy.


    def onclick(self, event):
        pass
       

    def onkey(self, event):
        if event.key == "d":
            print("removing last marker")
            self.lop = self.lop[:-1]
            handles = self.handles.pop(-1)
            handles.unrender()

        if event.inaxes and event.key == "w" and not self.in_area:
            self.lop += [(event.xdata, 0)]
            self.handles.append(handle_keeper())
            self.handles[-1].first_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="red")
            self.in_area = True
            plt.draw()

        elif event.inaxes and event.key == "w" and self.in_area:
            self.lop[-1] = (self.lop[-1][0], event.xdata)
            self.in_area = False
            self.handles[-1].second_line = plt.vlines(event.xdata, self.interval[0], self.interval[1], color="red")
            self.handles[-1].area = plt.fill_between(np.linspace(*self.lop[-1]), max(self.y), alpha=0.5, color="red")
            plt.draw()


        if event.key == "e":
            print("rerunning fit?")


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
