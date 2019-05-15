''' this takes a given atmosphere with structure in its oxygen abundance and
tries to recreate it with a fixed, well mixed oxygen abundance '''

'''
step 1: read in the given atmosphere
step 2: make an initial guess
for not converged:
    a) make a new pt file with new oxygen abundance
    b) run smart on new pt file
    c) calculate chi squared


'''
import smart
import numpy as np
import argparse
import matplotlib.pyplot as plt
import os
import scipy.optimize as sciopt
import shutil
import project_tools

'''
this script attempts to vary O2 abundance in a well mixed atmosphere to match
the O2 spectrum in an oxygen bulge scenario for a give band +/- 0.04 microns
'''


parser = argparse.ArgumentParser()
parser.add_argument('-experiment', dest='experiment_directory', type=str)
parser.add_argument('-init_o2', dest='init_o2', type=float)
parser.add_argument('-spec_type', dest='spec_type', type=str)
parser.add_argument('-band', dest='band', type=float)

args = parser.parse_args()


def rad_file(experiment_name, band):

    wlmin = band - 0.04
    wlmax = band + 0.04

    wnmax = 1e4 / wlmin
    wnmin = 1e4 / wlmax

    return '%s/case_%i_%icm_toa.rad' % (experiment_name, int(wnmin), int(wnmax))


def trnst_file(experiment_name, band):

    wlmin = band - 0.04
    wlmax = band + 0.04

    wnmax = 1e4 / wlmin
    wnmin = 1e4 / wlmax

    return '%s/case_%i_%icm.trnst' % (experiment_name, int(wnmin), int(wnmax))

def load_experiment_spectrum(experiment_name, band):

    if args.spec_type == 'transit':
        highres_fl = smart.readsmart.Trnst(trnst_file(experiment_name, band))
        # (r_planet/r_star)^2
#             highres_flux = highres_fl.tdepth

        # absolute radius
        highres_flux = highres_fl.absrad

    elif args.spec_type == 'direct':
        highres_fl = smart.readsmart.Rad(rad_file(experiment_name, band))

        highres_flux = highres_fl.pflux / highres_fl.sflux

    else:
        assert False

    lam = highres_fl.lam


    return highres_flux, lam

def chi_sq(experiment, model):

    return np.sum((experiment - model)**2)

def run_smart(ptfile, band, spec_type, place='degeneracy_test_output',
                R=100000., rm_when_done=False):
    NCPU=4
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

    sim.smartin.spec = '/Users/mcurr/RESEARCH/photochem_smart/fixed_input/specs/TRAPPIST-1.dat'
    sim.smartin.spec_skip = 5
    sim.smartin.spec_unit = 2 # Set unit to W/m^2/um

    sim.smartin.r_AU = 0.02916

    sim.smartin.radius = 0.915 * 6371.
    sim.smartin.grav = 0.930 * 9.8

    sim.lblin.radius = 0.915 * 6371.
    sim.lblin.grav = 0.930 * 9.8

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

    if spec_type == 'transit':
        # change to transit:
        sim.smartin.out_format = 9
        sim.smartin.itrnst = 2 # ray tracing
        sim.smartin.radstar = 0.121 # times sun radius
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

def model_spec(o2_abundance, experiment_name, band, spec_type, rm_when_done):

    pt_fl = project_tools.make_mixed_pt_fl(o2_abundance)


    project_tools.run_trappist(PT_DIR + pt_fl, band, wlrange=0.04, place='../data/smart_outputs/degeneracy/', R=100000., rm_when_done=rm_when_done)

    run_smart(PT_DIR+pt_fl, band,  spec_type, rm_when_done=rm_when_done)

    new_flux, new_lam = load_experiment_spectrum('../data/smart_outputs/degeneracy', band)
    shutil.rmtree('../data/smart_outputs/degeneracy')
    return new_flux, new_lam


def optimize_func(o2_abundance):
    new_flux, lam = model_spec(o2_abundance, 'degeneracy', args.band, args.spec_type, rm_when_done=False)

    return chi_sq(new_flux, experiment_flux)


PT_DIR = '../data/pt_fls/'

experiment_flux, lam = load_experiment_spectrum(args.experiment_directory, args.band)

res = sciopt.minimize(optimize_func, args.init_o2, bounds=[(0., 1.)])

print res
print res.x
print res.success
print res.message
