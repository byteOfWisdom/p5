import numpy as np
from matplotlib import pyplot as plt
from sys import argv
import std


class click_handler:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.interval = (min(y), max(y))
        self.lop = []


    def start(self):
        plt.connect("button_press_event", self.onclick)
        plt.connect("key_press_event", self.onkey)
        plt.plot(self.x, self.y)
        std.default.plt_pretty("Kanal", "Anzahl")
        plt.show()


    def onclick(self, event):
        pass

    def onkey(self, event):
        # if event.key == "d":
        #     print("removing last marker")
        #     self.lop = self.lop[:-1]

        if event.inaxes and event.key == "w":
            self.lop += [(event.xdata, event.ydata)]
            plt.vlines(event.xdata, self.interval[0], self.interval[1], color="red")
            plt.draw()

        if event.key == "e":
            print("rerunning fit?")


def let_user_click_peaks(x_values, y_values):
    handler = click_handler(x_values, y_values)
    handler.start()
    # ok might have wanted to do that earlier... let's not talk about that
    pass


def main():
    data = np.transpose(np.loadtxt(argv[1]))
    # out = argv[2]
    let_user_click_peaks(data[0], data[1])
    return None

if __name__ == "__main__":
    main()
