import project_tools
import numpy as np
o2_abundances = np.arange(0.05, 0.45, 0.05)
locs = ['upper', 'mixed']
for o2_abundance in o2_abundances:
    for loc in locs:
        project_tools.make_atmosphere(loc, o2_abundance)
