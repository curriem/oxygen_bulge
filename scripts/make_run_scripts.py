o2_bulge_options = ['upper', 'lower', 'mixed']

o2inv_options = ['0.21']
pt_shape_options = ['step', 'wedge']
n2_options = ['1.', '0.17', '0.007']
run_type_options = ['bands', 'full']

# upper o2 case

for o2_bulge in o2_bulge_options:

    with open('./job_%s.sh' % o2_bulge, 'wb') as f:
        for run_type in run_type_options:
            for o2inv in o2inv_options:
                for pt_shape in pt_shape_options:
                    for n2 in n2_options:
                        f.write('python ')
                        f.write('/Users/mcurr/RESEARCH/oxygen_bulge/scripts/')
                        f.write('run_experiment.py')
                        f.write(' -o2_loc %s ' % o2_bulge)
                        f.write(' -o2inv %s ' % o2inv)
                        f.write(' -n2scale %s ' % n2)
                        f.write(' -pt_shape %s ' % pt_shape)
                        f.write(' -run_type %s ' % run_type)
                        f.write('\n')
