import commands
import glob

systems = ['trappist', 'sun']
spec_types = ['transit', 'direct']
for system in systems:
    experiments = glob.glob('/Users/mcurr/RESEARCH/oxygen_bulge/data/%s/smart_outputs/*o2inv21*output' % system)
        for experiment in experiments:
            experiment = experiment.split('/')[-1]
            for spec in ['transit', 'direct']:
                print "python plot_report.py -experiment_name %s -system %s -spec_type %s" % (experiment, system, spec)
                commands.getoutput("python plot_report.py -experiment_name %s - system %s -spec_type %s" % (experiment, system, spec))
