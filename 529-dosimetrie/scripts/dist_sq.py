# %%
import std
import tomllib
from matplotlib import pyplot as plt
import numpy as np
import propeller as p


# %%
fhandle = open("../data/meta.toml", "rb")
meta = tomllib.load(fhandle)
fhandle.close()



def quadratic(x, a, c):
    denom = (x ** 2)
    return (a / denom) + c

# %%
fname = "../data/" + meta["abstand_sq"]["messung_1"]["csv_file"]
data = std.util.load_csv(fname, skiprows=1)

distance = p.ev(data[3] + 11.0, 0.5)
dose_before = p.ev(data[1], 0.1)
dose_after = p.ev(data[2], 0.1)
time = p.ev(data[0], 0.5)

dose_rate = (dose_after - dose_before) / time

print(r"\hline")
print(r"Abstand / cm & Dauer / s & $D_\text{1}$ / mSv & $D_\text{2}$ / mSv & Dosisleistung / mSv/s")

for i in range(len(dose_rate)):
    print(r"\num{", distance[i].format(), r"} & \num{", time[i].format(), r"} & \num{", dose_before[i].format(), r"} & \num{", dose_after[i].format(), r"} & \num{", dose_rate[i].format(), r"} \\")
print(r"\hline")


fit, (errs, rsq) = std.curve_fit(quadratic, distance, dose_rate)
params = p.ev(fit, errs)
print("\n parameter a & c:")
[print(params[i].format()) for i in range(len(fit))]

#plt.scatter(~distance, ~(1/np.sqrt(dose_rate)))
plt.errorbar(~distance, ~dose_rate, xerr = p.ve(distance)[1], yerr = p.ve(dose_rate)[1], label = "Messwerte", color = "tab:blue", **std.default.error_bar_def)
xrange = np.linspace(min(~distance) - 2, max(~distance) + 2)
plt.plot(xrange, quadratic(xrange, *fit), color = "tab:green", linewidth = 0.9, label = rf"Anpassungsfunktion, $R^2$={round(rsq,3)}")
plt.legend()
std.default.plt_pretty("Abstand / cm", "Dosisleistung / mSv/s")
plt.show()
