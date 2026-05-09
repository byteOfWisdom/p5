import std
from sys import argv 
import propeller as p


data = std.load_csv(argv[1], skiprows=1)
a = p.from_string("92.6664(19)")
b = p.from_string("154(19)")
print("A \t mu \t\t chi2")
for i in range(len(data)):
    print(f"{round(~data[0][i])} \t {((data[1][i] * a + b) / 1000).format()} keV \t {data[4][i]}")
