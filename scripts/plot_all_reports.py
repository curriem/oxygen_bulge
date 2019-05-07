import commands
import glob

experiments = glob.glob('/Users/mcurr/RESEARCH/oxygen_bulge/data/smart_outputs/*output')
for experiment in experiments:
    experiment = experiment.split('/')[-1]
    for spec in ['transit', 'direct']:
        print "python plot_report.py -experiment_name %s -spec_type %s" % (experiment, spec)
        commands.getoutput("python plot_report.py -experiment_name %s -spec_type %s" % (experiment, spec))
