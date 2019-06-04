import project_tools
import numpy as np
o2_abundance = 0.2
loc = 'upper'
cold_traps = np.logspace(-3, 0, 30)
for cold_trap in cold_traps:
    project_tools.make_atmosphere(loc, o2_abundance, cold_trap)
