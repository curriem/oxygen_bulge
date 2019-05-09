import project_tools

o2_bulge_options = ['upper', 'lower', 'mixed']

o2inv_options = ['0.21']
pt_shape_options = ['step', 'wedge']
n2_options = ['1.', '0.17', '0.007']
run_type_options = ['bands', 'full']

# upper o2 case
runfiles = []
for o2_bulge in o2_bulge_options:
    for run_type in run_type_options:
        for o2inv in o2inv_options:
            for pt_shape in pt_shape_options:
                for n2 in n2_options:
                    runfile = 'run_experiment.py -o2_loc %s -o2inv %s -n2scale %s -pt_shape %s -run_type %s' % (o2_bulge,
                                                                                                                       o2inv,
                                                                                                                       n2,
                                                                                                                       pt_shape,
                                                                                                                       run_type)
                    runfiles.append(runfile)

project_tools.write_slurm_script_python(runfiles,
                                        name='o2bulge',
                                        subname='o2bulge.sh',
                                        workdir='./',
					walltime='48:00:00')
