import smart
import os
import numpy as np
import argparse
import multiprocessing
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import project_tools


'''
how to run this script:

    1) choose an experiment name
        a) if you don't, a default name will be generated with specified run params, this is recommended
        b) usage: -exp <desired experiment name goes here>
    2) choose the location of your o2 bulge
        a) low -> bulge in lower atmosphere
        b) high -> bulge in upper atmosphere
        c) mixed -> mixed o2 in atmosphere
        d) usage: -o2_loc <low, high, or mixed (choose ONE)>
    3) choose whether you want o2 o2 collisionally induced absorption on or off
        a) use -o4cia_on for on
        b) use -o4cia_off for off
    4) choose whether you want water on or off
        a) use -h2o_on for on
        b) use -h2o_off for off
    5) choose whether you want transit or direct spectroscopy
        a) use -transit for transit
        b) use -direct for direct imaging
    6) set maximum o2 abundance
        a) usage -o2inv <mixing ratio between 0 and 1 here>
        b) examples:
            i) for a lower atmosphere o2 bulge case, this will be the o2 abundance at the surface
            ii) for an upper atmosphere o2 bulge case, this will be the o2 abundance at the highest point in the atmosphere
            iii) for a mixed atmosphere, this will be the o2 abundance at all pressures
    7) set the resolution you wish to run with
        a) usage: -R <resolution you want to run with>

    EXAMPLE RUNS:
    python run_experiment.py -o2_loc low -o4cia_on -h2o_on -transit -o2inv 0.21 -R 100000
    python run_experiment.py -o2_loc high -o4cia_on -h2o_off -direct -o2inv 0.4 -R 150
    python run_experiment.py -o2_loc mixed -o4cia_off -h2o_off -direct -o2inv 0.01 -R 15000
    python run_experiment.py -exp experiment1 -o2_loc high -o4cia_off -h2o_off -direct -o2inv 0.21 -R 10000

'''

parser = argparse.ArgumentParser()
#parser.add_argument('-exp', dest='experiment_name', default=None, type=str, help='Name of experiment')
parser.add_argument('-o2_loc', dest='o2_loc', type=str, help='upper, lower, mixed')
parser.add_argument('-trop_loc', dest='trop_loc', default=None, type=float, help='pressure [bar] which separates upper and lower atmosphere regimes')
#parser.add_argument('-o4cia_on', dest='o2o2', action='store_true', help='turn on o4 cia')
parser.add_argument('-o4cia_off', dest='o2o2', action='store_false', help='turn off o4 cia')
#parser.add_argument('-h2o_on', dest='h2o', action='store_true', help='turn on water')
parser.add_argument('-h2o_off', dest='h2o', action='store_false', help='turn off water')
parser.add_argument('-o2inv', dest='o2inv', default=0.21, type=float, help='o2 abundance in upper, lower, and mixed cases')
parser.add_argument('-n2scale', dest='n2_scalar', default=1., type=float, help='times PAL N2')
parser.add_argument('-R', dest='resolution', default=100000, type=int, help='Resolution of spectra')
parser.add_argument('-pt_shape', dest='pt_shape', default='wedge', type=str)
parser.add_argument('-run_type', dest='run_type', type=str)

args = parser.parse_args()

# default things:

experiment_name = args.o2_loc

if not args.h2o:
    h2o = False
    experiment_name += '_H2Ooff'
else:
    h2o = True

if not args.o2o2:
    o2o2 = False
    experiment_name += '_O4off'
else:
    o2o2 = True

if args.trop_loc == None:
    trop_loc = 0.1
else:
    trop_loc = args.trop_loc
    experiment_name = experiment_name + '_trop%i'%int(trop_loc*10)

experiment_name = experiment_name + '_o2inv%s' % str(int(100*args.o2inv)).zfill(2)

if args.n2_scalar == 1.:
    experiment_name = experiment_name + '_n2scale1'
else:
    experiment_name = experiment_name + '_n2scale%s' % str(args.n2_scalar).split('.')[-1]

experiment_name = experiment_name + '_%s' % args.pt_shape

experiment_name = experiment_name + '_%s' % args.run_type


PT_DIR = '../data/pt_fls/'
SMART_OUTPUTS = '../data/smart_outputs/'

NCPU = multiprocessing..cpu_count() - 2

# experiment_name = 'o2%s%i_trop%i_water%i_o4cia%i_%s' % (args.o2_loc,
#                                                         int(100*args.o2_abundance),
#                                                         int(args.trop_loc*10),
#                                                         int(args.h2o),
#                                                         int(args.o2o2),
#                                                         args.pt_shape)
# else:
#     experiment_name = args.experiment_name

print '\n'
print 'Experiment name:', experiment_name
# print 'Location of O2 bulge:', args.o2_loc
# print 'Location of tropopause:', args.trop_loc
# print 'O2 O2 CIA on?', args.o2o2
# print 'Water on?', args.h2o
# print 'Transit spectroscopy?', args.transit
# print 'Maximum O2 abundance:', args.o2_abundance
# print 'Resolution:', args.resolution
# print 'PT shape:', args.pt_shape
print '\n'


OUTPUT_DIR = SMART_OUTPUTS + experiment_name + '_output'


# check if pt file exists for this experiment
pt_fl = experiment_name + '.pt'
#pt_fl = 'o2%s%i_water%i.pt' % (args.o2_loc, int(100*args.o2_abundance), int(args.h2o))
if not os.path.isfile(PT_DIR + pt_fl):
    # make pt file if it does not exist already
    pt_fl = project_tools.make_pt_fl(args.o2inv, args.o2_loc, trop_loc, h2o,
                                args.pt_shape, experiment_name,
                                N2_scalar=args.n2_scalar)


if args.run_type == 'bands':
    bands_to_run = [0.63, 0.68, 0.76, 1.27]
    wlrange = 0.04
elif args.run_type == 'full':
    bands_to_run = [0.95]
    wlrange = 0.55

for band in bands_to_run:
    #print "run_smart({PT_DIR} + pt_fl, band, place=OUTPUT_DIR, R=args.resolution, o2o2=args.o2o2, transit=args.transit)"
    project_tools.run_trappist(PT_DIR + pt_fl, band=band, wlrange=wlrange,
                            place=SMART_OUTPUTS + '/trappist/' + experiment_name + '_output', R=args.resolution,
              o2o2=o2o2, NCPU=NCPU, addn2=False)
    project_tools.run_sun(PT_DIR + pt_fl, band=band, wlrange=wlrange,
                            place=SMART_OUTPUTS + '/sun/' + experiment_name + '_output', R=args.resolution,
              o2o2=o2o2, NCPU=NCPU, addn2=False)
