#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 10:08:52 2019

@author: mcurr
"""

import matplotlib.pyplot as plt
import matplotlib
import argparse
import coronagraph as cg
import smart
import numpy as np

############ add arguments #####################

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
parser.add_argument('-highres_R', dest='highres_R', default=100000, type=int, help='Resolution of spectra')
parser.add_argument('-lowres_R', dest='lowres_R', default=150, type=int, help='Resolution of spectra')


args = parser.parse_args()

##############################################

####### fix plotting default values #######

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 18}

figure = {'figsize' : (12, 10)}

lines = {'linewidth' : 3}

matplotlib.rc('font', **font)
matplotlib.rc('figure', **figure)
matplotlib.rc('lines', **lines)

#############################################

################ plotting functions #################3

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
    
def rad_file(experiment_name, band):
    
    wlmin = band - 0.04
    wlmax = band + 0.04
    
    wnmax = 1e4 / wlmin
    wnmin = 1e4 / wlmax
    
    return '../data/smart_outputs/%s_output/case_hitran2012_%i_%icm_toa.rad' % (experiment_name, int(wnmin), int(wnmax))


def trnst_file(experiment_name, band):
    
    wlmin = band - 0.04
    wlmax = band + 0.04
    
    wnmax = 1e4 / wlmin
    wnmin = 1e4 / wlmax
    
    return '../data/smart_outputs/%s_output/case_hitran2012_%i_%icm.trnst' % (experiment_name, int(wnmin), int(wnmax))

def degrade_spectrum(flux, lam, band, R, span=0.04):
    lam_min = band - span
    lam_max = band + span
    # Construct new low-res wavelength grid
    wl, dwl = cg.noise_routines.construct_lam(lam_min, lam_max, Res=R)
    flr = cg.downbin_spec(flux, lam, wl, dlam=dwl)
    return wl, flr

def plot_spectra(experiment_name, band, highres_wl, highres_flux,
                 lowres_wl, lowres_flux,
                 highres_R, lowres_R, transit=False):
    plt.figure()
    plt.title(r'%s $\mu$m band' % str(band))
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    if transit:
        #plt.ylabel(r'Transit Depth [(R$_{\mathrm{planet}}$ / R$_{\mathrm{star}}$)$^2$]')
        plt.ylabel('Planet Radius [km]')
    else:
        #plt.ylabel(r'Flux [W/m$^2$/$\mu$m]')
        plt.ylabel(r'Reflectivity')
    
    plt.plot(highres_wl, highres_flux, color='darkslategray', label='mixed, R=%i' % highres_R)
    plt.plot(lowres_wl, lowres_flux, color='paleturquoise', linestyle='dashed', label='mixed, R=%i' % lowres_R)
    
    plt.legend()
    #plt.show()
    plt.savefig('../plots/%s_%s.pdf' % (experiment_name, str(int(band*100)).zfill(3)))
    
def extract_and_plot(experiment_name):
    for band in [0.63, 0.68, 0.76, 1.27]:

        if args.transit:
            highres_fl = smart.readsmart.Trnst(trnst_file(experiment_name, band))
            
            
            # (r_planet/r_star)^2
#             highres_flux = highres_fl.tdepth
            
            # absolute radius
            highres_flux = highres_fl.absrad

            
        else:
            highres_fl = smart.readsmart.Rad(rad_file(experiment_name, band))
            
            highres_flux = highres_fl.pflux / highres_fl.sflux
            
        highres_wl = highres_fl.lam
        lowres_wl, lowres_flux = degrade_spectrum(highres_flux, highres_wl,
                                                   band, R=args.lowres_R)

        plot_spectra(experiment_name, band, highres_wl, highres_flux,
                     lowres_wl, lowres_flux,
                     args.highres_R, args.lowres_R, transit=args.transit)


#####################################################

############# get experiment name ##############

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
    
################################################
    
############# point to directories ##############
    
PT_DIR = '../data/pt_fls'
SMART_OUTPUTS = '../data/smart_outputs'
OUTPUT_DIR = SMART_OUTPUTS + experiment_name + '_output'

##################################################

extract_and_plot(experiment_name)
    

    

