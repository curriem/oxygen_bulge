import project_tools
import numpy as np
import glob
import sys

name = sys.argv[1]
star = sys.argv[2]
planet = sys.argv[3]

pt_files1 = glob.glob('../data/pt_fls/*_[1,2,3,4,5]*')
pt_files2 = glob.glob('../data/pt_fls/*_01*')
pt_files3 = glob.glob('../data/pt_fls/*_02*')
pt_files4 = glob.glob('../data/pt_fls/*_03*')
pt_files5 = glob.glob('../data/pt_fls/*_04*')
pt_files6 = glob.glob('../data/pt_fls/*_05*')
pt_files7 = glob.glob('../data/pt_fls/*_06*')
pt_files8 = glob.glob('../data/pt_fls/*_07*')
pt_files9 = glob.glob('../data/pt_fls/*_08*')

new = glob.glob('../data/pt_fls/*_09*') +  glob.glob('../data/pt_fls/*_1*')

pt_files_master = [pt_files9, pt_files8, pt_files7, pt_files6]#, pt_files5, pt_files6]

pt_fls = glob.glob('../data/pt_fls/*_00*')
pt_files_master = [pt_fls]

for n in range(1):
        pt_files = pt_files_master[n]
	runfiles = []
	for pt_file in pt_files:
	    pt_file = pt_file.split('/')[-1]
	    print pt_file
	    runfile = 'run_experiment.py -star %s -planet %s -pt_fl %s -wlmin 0.04 -wlmax 1.5' % (star, planet, pt_file)
	    runfiles.append(runfile)


	project_tools.write_slurm_script_python(runfiles,
						name=name+str(n+1),
						subname=name+str(n+1)+'.sh',
						workdir='./',
						walltime='100:00:00')
