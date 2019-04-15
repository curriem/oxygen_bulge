import smart
import os
import numpy as np
import argparse
import multiprocessing
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


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
parser.add_argument('-exp', dest='experiment_name', default=None, type=str, help='Name of experiment')
parser.add_argument('-o2_loc', dest='o2_loc', type=str, help='high, low, mixed')
parser.add_argument('-o4cia_on', dest='o2o2', action='store_true', help='turn on o4 cia')
parser.add_argument('-o4cia_off', dest='o2o2', action='store_false', help='turn off o4 cia')
parser.add_argument('-h2o_on', dest='h2o', action='store_true', help='turn on water')
parser.add_argument('-h2o_off', dest='h2o', action='store_false', help='turn off water')
parser.add_argument('-transit', dest='transit', action='store_true', help='use for direct imaging spectroscopy')
parser.add_argument('-direct', dest='transit', action='store_false', help='use for transit spectroscopy')
parser.add_argument('-o2inv', '--oxygenabundance', dest='o2_abundance', default=0.21, type=float, help='o2 abundance in upper, lower, and mixed cases')
parser.add_argument('-R', '-resolution', dest='resolution', default=100000, type=int, help='Resolution of spectra')

args = parser.parse_args()



def run_smart(ptfile, band, place='output', R=100000., rm_when_done=False,
              o2o2=True, transit=False):
    '''    
    Parameters
    -----------
    place : str
        Directory to create to run this test in
    rm_when_done : bool
        Set to remove ``place'' when the simulation finishes
    R : float
        resolving power
    '''
    try:
        os.mkdir(place)
    except OSError:
        pass
    
    # load pt file (relative import frim photochem directory)
    atm = smart.interface.Atmosphere.from_pt(ptfile, atm_dir=place)
    
    # Set atmosphere parameters
    atm.set_dataframe_from_profiles()
    atm.write_atm()
    
    # Prep simulation
    sim = smart.interface.Smart(tag='case', atmosphere=atm)
    sim.smartin.source = 3 # solar and thermal sources
    sim.smartin.out_format = 1 # no transit calculations
    
    sim.set_run_in_place(place=place)
    
    wlmin = band - 0.04
    wlmax = band + 0.04
#     wlmin = band - 0.1
#     wlmax = band + 0.1
    
    # change the wavelength range
    sim.smartin.minwn = 1e4 / wlmax
    sim.smartin.maxwn = 1e4 / wlmin
    sim.lblin.minwn = 1e4 / wlmax
    sim.lblin.maxwn = 1e4 / wlmin
    
    # watch out for indexing: check if name equals the gas I want
    if not o2o2:
        sim.atmosphere.gases[6].cia_file = None
    
    if transit: 
        # change to transit:
        sim.smartin.out_format = 9
        sim.smartin.itrnst = 2 # ray tracing 
        sim.smartin.radstar = 1 # times sun radius
        sim.smartin.irefract = 1 # turn on refraction     
    
    # modify resolution of spectrum 
    delt_nu = 1e4 / ((wlmax - wlmin)/2 * R)
    sim.smartin.FWHM = delt_nu
    sim.smartin.sample_res = delt_nu
    
    # set to run on this machine
    sim.set_executables_automatically()
    
    # generate LBLABC scripts
    sim.gen_lblscripts()
    
    # run LBLABC on all but two CPUs on machine
    sim.run_lblabc(ncpu=NCPU)
    
    #sim.tag = sim.tag + '_' + tag + '_' + 'R' + str(R)
    
    # run SMART
    sim.write_smart(write_file=True)
    sim.run_smart()
    
    # remove temp directory
    if rm_when_done:
        os.system('rm -r %s' % place)
        
    return

def make_pt_fl(max_o2_abundance, o2_location, h2o, experiment_name):
    # Standard Earth O_2 mixing ratio:
    earth_standard_ptfile = PT_DIR + '/earth_standard_icrccm_vmix.pt'
    header = "Press        Temp       H2O             CO2        O3            N2O          CO           CH4         O2"
    earth_pt_data = np.genfromtxt(earth_standard_ptfile, skip_header=1)
    earth_pt_data_inv = earth_pt_data.T
    earth_pressure = earth_pt_data_inv[0, :]
    earth_o2_mixing_ratio = earth_pt_data_inv[8, :]
    
    #new_pt_fl_name = 'o2%s%s_water%i.pt' % (o2_location, str(int(max_o2_abundance*100)).zfill(2), int(h2o))
    new_pt_fl_name = experiment_name + '.pt'
    pressure = [0, -1, -7]
    if o2_location == 'low':
        o2_mix = [np.log10(max_o2_abundance), -5, -14]
        fill_val = np.log10(max_o2_abundance)
    elif o2_location == 'high':
        o2_mix = [-14, -5, np.log10(max_o2_abundance)]
        fill_val = -14
    else:
        # mixed o2 abundance
        o2_mix = [np.log10(max_o2_abundance)]*3
        fill_val = np.log10(max_o2_abundance)
    
    func = interp1d(pressure, o2_mix, bounds_error=False,
                    fill_value=fill_val)
    new_o2_mix = func(np.log10(earth_pressure))
    new_o2_mix = 10**new_o2_mix
    
    new_pt_inv = np.copy(earth_pt_data_inv)
    new_pt_inv[8, :] = new_o2_mix
    
    if not h2o:
        new_pt_inv[2, :] /= 1000
        
    new_pt_data = new_pt_inv.T
    np.savetxt(PT_DIR+new_pt_fl_name, new_pt_data, delimiter='   ', header=header, comments='')
    
    return new_pt_fl_name
    
    
def plot_pt(pt_fl):
    pt = np.genfromtxt(PT_DIR + pt_fl, skip_header=1)
    pt_data = pt.T
    pressure = pt_data[0, :]
    temperature = pt_data[1, :]
    o2 = pt_data[8, :]
    plt.figure()
    plt.loglog(o2, pressure)
    plt.xlabel(r'O$_2$ mixing ratio')
    plt.ylabel('Pressure [bar]')
    plt.show()
    

    
if __name__ == '__main__':
    PT_DIR = '../data/pt_fls/'
    SMART_OUTPUTS = '../data/smart_outputs/'

    NCPU = 1

    print 'Experiment name:', args.experiment_name
    print 'Location of O2 bulge:', args.o2_loc
    print 'O2 O2 CIA on?', args.o2o2
    print 'Water on?', args.h2o
    print 'Transit spectroscopy?', args.transit
    print 'Maximum O2 abundance:', args.o2_abundance
    print 'Resolution:', args.resolution
    print '\n'
    
    if args.experiment_name == None:
        if args.transit:
            spec = 'transit'
        else:
            spec = 'direct'
        experiment_name = 'o2%s%i_water%i_o4cia%i_%s' % (args.o2_loc,
                                                         int(100*args.o2_abundance),
                                                         int(args.h2o),
                                                         int(args.o2o2),
                                                         spec)
    else:
        experiment_name = args.experiment_name
    
    OUTPUT_DIR = SMART_OUTPUTS + experiment_name + '_output'
        
    
    
    # check if pt file exists for this experiment
    pt_fl = experiment_name + '.pt'
    #pt_fl = 'o2%s%i_water%i.pt' % (args.o2_loc, int(100*args.o2_abundance), int(args.h2o))
    if not os.path.isfile(PT_DIR + pt_fl):
        # make pt file if it does not exist already
        make_pt_fl(args.o2_abundance, args.o2_loc, args.h2o, experiment_name)

    
    bands_to_run = [0.63, 0.68, 0.76, 1.27]
    
    for band in bands_to_run:
        #print "run_smart({PT_DIR} + pt_fl, band, place=OUTPUT_DIR, R=args.resolution, o2o2=args.o2o2, transit=args.transit)"
        run_smart(PT_DIR + pt_fl, band, place=OUTPUT_DIR, R=args.resolution, 
                  o2o2=args.o2o2, transit=args.transit)
    
    
