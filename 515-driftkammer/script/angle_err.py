import std
import numpy as np
import propeller as p

a = p.ev(108, 5) 
s0 = p.ev(22.5 * 8.5, 25)
n = np.arange(0, 50, 2)
d = n * 8.5 - s0

theta = 180. * np.atan(d / a) / np.pi

for i in range(len(theta)):
    print(f"${n[i]}$ & $\\num" + "{" +  theta[i].format() + "}$ \\\\")
