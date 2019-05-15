import project_tools
import numpy as np
import glob
import sys

name = sys.argv[1]

pt_files = glob.glob('../data/pt_fls/*')
locs = ['upper', 'lower']

star = 'trappist1'
planet = 't1e'

runfiles = []
for pt_file in pt_files:
    runfile = 'run_experiment.py -star %s -planet %s -pt_fl %s -wlmin 0.04 -wlmax 1.5' % (star, planet, pt_file)
    runfiles.append(runfile)


project_tools.write_slurm_script_python(runfiles,
                                        name=name,
                                        subname=name+'.sh',
                                        workdir='./',
					walltime='100:00:00')
