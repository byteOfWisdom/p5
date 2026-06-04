import std
import numpy as np
import propeller as p

muon_file = "../data/rate_myonen.csv"
electron_file = "../data/rate_electronen.csv"

muon_data = std.load_csv(muon_file, ",", 1)
electron_data = std.load_csv(electron_file, ",", 1)

ekey = np.argsort(electron_data[2])
e_counts = electron_data[0][ekey]
e_times = electron_data[1][ekey]
e_times_w_err = p.ev(e_times,e_times*3e-2)

mkey = [np.where(muon_data[2] == e_voltage)[0][0] for e_voltage in electron_data[2][ekey]]
m_counts = muon_data[0][mkey]
m_times = muon_data[1][mkey]
m_times_w_err = p.ev(m_times, m_times*3e-2)

#for r, u in zip(muon_data[0] / muon_data[1], muon_data[2]):
    #print(u, r)

voltage = electron_data[2][ekey]
voltage_w_err = p.ev(voltage,voltage*1e-2)

muon_rate = m_counts / m_times_w_err
electron_rate = e_counts / e_times_w_err
ratio = electron_rate / muon_rate

print(r"$U_\text{Szint}$ / kV & $n_\mu$ / [o.E.] & $t_\mu$ / s & $n_e$ / [o.E.] & $t_e$  / s \\")
print(r"\hline")
for i in range(len(voltage_w_err)):
    print(r"$\num{", voltage_w_err[i].format(), r"}$&", r"$\num{", m_counts[i], r"}$&", r"$\num{", m_times_w_err[i].format(), r"}$&", r"$\num{", e_counts[i], r"}$&", r"$\num{", e_times_w_err[i].format(), r"}$ \\")
print(r"\hline")

# print(r"$U_\text{Szint}$ / kV & $a_\mu$ / $\unit{\per\s}$ & $a_e$ / $\unit{\per\s}$ & A [o.E.] \\")
# print(r"\hline")
# for i in range(len(voltage_w_err)):
#     print(r"$\num{", voltage_w_err[i].format(), r"}$&",r"$\num{", muon_rate[i].format(), r"}$&", r"$\num{", electron_rate[i].format(), r"}$&", r"$\num{", ratio[i].format(),r"}$ \\")
# print(r"\hline")

# for i in range(len(e_counts)):
#     print("e count:", e_counts[i])
#     print("e time:", e_times[i])

# std.print_csv_table({
#                         "voltage": voltage,
#                         "muon rate": muon_rate,
#                         "electron rate": electron_rate,
#                         "ratio": ratio
#                     }, "../figs/ratios.csv")
