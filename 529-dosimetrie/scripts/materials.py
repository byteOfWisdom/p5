import std
# import propeller as p
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

absorb_coeff_a = - np.log(transmission_a)
absorb_coeff_b = - np.log(transmission_b)

std.default.plt_errorbar(material_atomic_nums[1:], absorb_coeff_a[1:], "Ohne Filter")
std.default.plt_errorbar(material_atomic_nums[1:], absorb_coeff_b[1:], "Mit Filter")
std.default.plt_finish("Ordnungszahl", "$\\mu$")
