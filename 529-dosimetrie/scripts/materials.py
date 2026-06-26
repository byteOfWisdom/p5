import std
import propeller as p
import numpy as np

import shielding

params_a = shielding.metadata["messung_3"]
params_b = shielding.metadata["messung_4"]

print(params_a)
print(params_b)

transmission_a = shielding.load_and_normalize(params_a["csv_file"])
transmission_b = shielding.load_and_normalize(params_b["csv_file"])

materials = ["Leer", "C", "Al", "Fe", "Cu", "Zr", "Ag"]
material_atomic_nums = [0, 6, 13, 26, 29, 40, 47]
literature_coeffs = [0, 0, 0.929, 0, 0, 0, 0]

mu_a = p.ev(0.850, 0.046)
mu_b = p.ev(0.828, 0.024)

thickness_measured_a = transmission_a[2] / mu_a
thickness_measured_b = transmission_a[2] / mu_b
thickness_true = transmission_a[2] / 0.929
print("with measured mu: ", thickness_measured_a.format(), "mm")
print("with measured mu: ", thickness_measured_b.format(), "mm")
print("with true mu: ", thickness_true.format(), "mm")

absorb_coeff_a = - np.log(transmission_a) / thickness_measured_a
absorb_coeff_b = - np.log(transmission_b) / thickness_measured_a

std.default.plt_errorbar(material_atomic_nums[1:], absorb_coeff_a[1:], "Ohne Filter")
std.default.plt_errorbar(material_atomic_nums[1:], absorb_coeff_b[1:], "Mit Filter")
std.default.plt_errorbar(material_atomic_nums[1:], literature_coeffs[1:], "Literatur")
std.default.plt_finish("Ordnungszahl", "$\\mu$ / $mm^{-1}$")
