import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

MAT = scipy.io.loadmat('../data/abiotic_O2_2014/results.mat')
T = [MAT['T'][:, 0], MAT['T'][:, 1], MAT['T'][:, 2]]
p = [MAT['p'][:, 0], MAT['p'][:, 1], MAT['p'][:, 2]]

fwater = [MAT['f'][:, 0], MAT['f'][:, 1], MAT['f'][:, 2]]

inds_to_remove = [2, 7, 12, 17,
                  22, 27, 32, 37,
                  42, 47, 52, 57,
                  62, 67, 72, 77]
bool_inds_to_remove = np.ones(80, dtype=bool)
bool_inds_to_remove[inds_to_remove] = False


n2pal = [1, 0.17, 0.007] # times earth PAL

fco2 = [1e-3, 0.1, 0.9] # mol / mol


header = ["Press",
          "Temp",
          "H2O",
          "CO2",
          "O2",
          "N2"]

header = "        ".join(header)

earth_pt = '../data/pt_fls/earth_standard_icrccm_vmix.pt'
earth_pt_data = np.genfromtxt(earth_pt, skip_header=1)
earth_sum_gas_fractions = np.sum(earth_pt_data[:, 2:], axis=1)
earth_n2_fraction = np.ones_like(earth_sum_gas_fractions) - earth_sum_gas_fractions
earth_pressure = earth_pt_data[:, 0]
earth_n2_pressure = earth_pressure * earth_n2_fraction



func = interp1d(earth_pressure, earth_n2_pressure)

n2_press_case = func(p[0])


n2_frac_case = n2_press_case / p[0]


earth_n2_fraction =  (1 - fwater[0][0] - fco2[0]) * np.ones(80)
earth_n2_pressure = earth_n2_fraction * p[0]

for n in range(3):
    fn2 = n2pal[n] * earth_n2_pressure / p[n]
    fn2 = fn2[bool_inds_to_remove]
    T[n] = T[n][bool_inds_to_remove]
    p[n] = p[n][bool_inds_to_remove]
    fwater[n] = fwater[n][bool_inds_to_remove]


    print fn2 + fwater[n] + fco2[n]*np.ones(64)

    fo2 = np.ones_like(fn2) - fn2 - fwater[n] - fco2[n]
