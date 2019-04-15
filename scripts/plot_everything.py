import commands

spec_options = ['-transit', '-direct']
water_options = ['-h2o_on', '-h2o_off']
bulge_options = ['high', 'low', 'mixed']
o2inv_options = ['0.21', '0.40']

for spec in spec_options:
    for water in water_options:
        for bulge in bulge_options:
            for o2_inv in o2inv_options:
                commands.getoutput('python plot_experiment.py -o2_loc %s -o4cia_on %s %s -o2inv %s -highres_R 100000 -lowres_R 150' % (bulge, water, spec, o2_inv))

