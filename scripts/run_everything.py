import commands

spec_options = ['-transit', '-direct']
water_options = ['-h2o_on', '-h2o_off']
bulge_options = ['high', 'low', 'mixed']
#o2inv_options = ['0.1', '0.2', '0.3', '0.4']
o2inv_options = ['0.00001']
for spec in spec_options:
    for water in water_options:
        for bulge in bulge_options:
            for o2_inv in o2inv_options:
                commands.getoutput('python run_experiment.py -o2_loc %s -o4cia_on %s %s -o2inv %s -R 100000' % (bulge, water, spec, o2_inv))

