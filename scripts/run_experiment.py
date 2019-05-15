import smart
import os
import numpy as np
import argparse
import multiprocessing
import matplotlib.pyplot as plt
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
parser.add_argument('-planet', dest='planet', type=str)
parser.add_argument('-star', dest='star', type=str)
parser.add_argument('-pt_fl', dest='pt_fl', type=str)
parser.add_argument('-wlmin', dest='wlmin', type=float)
parser.add_argument('-wlmax', dest='wlmax', type=float)
parser.add_argument('-R', dest='resolution', default=100000, type=int, help='Resolution of spectra')


#parser.add_argument('-exp', dest='experiment_name', default=None, type=str, help='Name of experiment')
# parser.add_argument('-o2_loc', dest='o2_loc', type=str, help='upper, lower, mixed')
# parser.add_argument('-trop_loc', dest='trop_loc', default=None, type=float, help='pressure [bar] which separates upper and lower atmosphere regimes')
# #parser.add_argument('-o4cia_on', dest='o2o2', action='store_true', help='turn on o4 cia')
# parser.add_argument('-o4cia_off', dest='o2o2', action='store_false', help='turn off o4 cia')
# #parser.add_argument('-h2o_on', dest='h2o', action='store_true', help='turn on water')
# parser.add_argument('-h2o_off', dest='h2o', action='store_false', help='turn off water')
# parser.add_argument('-o2inv', dest='o2inv', default=0.21, type=float, help='o2 abundance in upper, lower, and mixed cases')
# parser.add_argument('-n2scale', dest='n2_scalar', default=1., type=float, help='times PAL N2')
# parser.add_argument('-R', dest='resolution', default=100000, type=int, help='Resolution of spectra')
# parser.add_argument('-pt_shape', dest='pt_shape', default='wedge', type=str)
# parser.add_argument('-run_type', dest='run_type', type=str)
# parser.add_argument('-pt_only', dest='pt_only', action='store_false')

args = parser.parse_args()




PT_DIR = '../data/pt_fls/'
SMART_OUTPUTS = '../data/smart_outputs/'


project_tools.run_smart(args.star, args.planet, args.pt_fl,
                        wlmin=args.wlmin, wlmax=args.wlmax, NCPU=1)
