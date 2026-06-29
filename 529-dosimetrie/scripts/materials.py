import std
import propeller as p
import numpy as np
import xraydb

import shielding

params_a = shielding.metadata["messung_3"]
params_b = shielding.metadata["messung_4"]

print(params_a)
print(params_b)

transmission_a = shielding.load_and_normalize(params_a["csv_file"])
transmission_b = shielding.load_and_normalize(params_b["csv_file"])

materials = ["Leer", "C", "Al", "Fe", "Cu", "Zr", "Ag"]
material_atomic_nums = [0, 6, 13, 26, 29, 40, 47]

mu_a = p.ev(0.850, 0.046)
mu_b = p.ev(0.828, 0.024)

thickness_measured_a = - np.log(transmission_a[2]) / mu_a
thickness_measured_b = - np.log(transmission_a[2]) / mu_b
thickness_true = transmission_a[2] / (0.1 * xraydb.material_mu("Al", 20e3))

def ev_mu(elem):
    value = 0.1 * xraydb.material_mu(elem, 20e3)
    edges = 0.1 * np.array([xraydb.material_mu(elem, e) for e in [15e3, 25e3]])
    # value = 0.1 * xraydb.material_mu(elem, 17.3e3)
    # edges = 0.1 * np.array([xraydb.material_mu(elem, e) for e in [17.1e3, 17.5e3]])
    diff = np.abs(edges - value)
    return p.ev(value, np.average(diff))


literature_coeffs = [0] + [ev_mu(elem) for elem in materials[1:]]

print(literature_coeffs[-2])
print("with measured mu: ", thickness_measured_a.format(), "mm")
print("with measured mu: ", thickness_measured_b.format(), "mm")
print("with true mu: ", thickness_true.format(), "mm")

absorb_coeff_a = - np.log(transmission_a) / thickness_measured_a
absorb_coeff_b = - np.log(transmission_b) / thickness_measured_a

# absorb_coeff_a = - np.log(transmission_a) / thickness_true
# absorb_coeff_b = - np.log(transmission_b) / thickness_true

std.default.plt_errorbar(material_atomic_nums[1:], absorb_coeff_a[1:], "Ohne Filter", "D")
std.default.plt_errorbar(material_atomic_nums[1:], absorb_coeff_b[1:], "Mit Filter", "D")
std.default.plt_errorbar(material_atomic_nums[1:], literature_coeffs[1:], "Literatur", "D")
std.default.plt_finish("Ordnungszahl", "$\\mu$ / $mm^{-1}$")

ids = np.arange(len(materials) - 1)
std.default.plt_errorbar(ids + 0.1, absorb_coeff_a[1:], "Ohne Filter", "D")
std.default.plt_errorbar(ids - 0.1, absorb_coeff_b[1:], "Mit Filter", "D")
std.default.plt_errorbar(ids, literature_coeffs[1:], "Literatur", "D")
std.default.plt.gca().set_xticks(ids, materials[1:])
std.default.plt.yscale("log")
std.default.plt_finish("Element", "$\\mu$ / $mm^{-1}$")


table = {
    "Element": materials,
    "Ordnungszahl": material_atomic_nums,
    "$\\mu_\\text{ohne Filter}$ / $mm^{-1}$": absorb_coeff_a,
    "$\\mu_\\text{mit Filter}$ / $mm^{-1}$": absorb_coeff_b,
}

std.print_tex_table(table, "../latex/material_mu.table")
