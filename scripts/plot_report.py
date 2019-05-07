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
parser.add_argument('-experiment_name', dest='experiment', type=str)
parser.add_argument('-spec_type', dest='spec_type', type=str)
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
args = parser.parse_args()

lowres_R = 150
highres_R = 100000

PT_DIR = '../data/pt_fls/'
SMART_OUTPUTS = '../data/smart_outputs/'



experiment = args.experiment



experiment_split = experiment.split('_')

experiment = '_'.join(experiment_split[:-1])

if len(experiment_split) == 6:
    o2_loc = experiment_split[0]
    o2_inv = float(experiment_split[1].split('v')[-1]) / 100
    n2scale = experiment_split[2]
    pt_shape = experiment_split[3]
    run_type = experiment_split[4]

####### fix plotting default values #######

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

figure = {'figsize' : (18, 10)}

lines = {'linewidth' : 1}

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

def rad_file(experiment_name, band, wlspan=0.04):

    wlmin = band - wlspan
    wlmax = band + wlspan

    wnmax = 1e4 / wlmin
    wnmin = 1e4 / wlmax

    return '../data/smart_outputs/%s_output/case_hitran2012_%i_%icm_toa.rad' % (experiment_name, int(wnmin), int(wnmax))


def trnst_file(experiment_name, band, wlspan):

    wlmin = band - wlspan
    wlmax = band + wlspan

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

def extract(experiment_name, band):

    if args.spec_type == 'transit':
        highres_fl = smart.readsmart.Trnst(trnst_file(experiment_name, band, 0.04))


        # (r_planet/r_star)^2
#             highres_flux = highres_fl.tdepth

        # absolute radius
        highres_flux = highres_fl.absrad


    elif args.spec_type == 'direct':
        highres_fl = smart.readsmart.Rad(rad_file(experiment_name, band, 0.04))

        highres_flux = highres_fl.pflux / highres_fl.sflux

    highres_wl = highres_fl.lam
    lowres_wl, lowres_flux = degrade_spectrum(highres_flux, highres_wl,
                                               band, R=lowres_R)

    return highres_wl, highres_flux, lowres_wl, lowres_flux

def get_pt_data(experiment_name):

    pt_data = np.genfromtxt(PT_DIR + experiment_name + '.pt', skip_header=1)
    pressure = pt_data.T[0, :]
    abundance = pt_data.T[8, :]

    return pressure, abundance

def extract_pt_dict(pt_fl):

    header = "P        T       H2O             CO2        O3            N2O          CO           CH4         O2     N2"
    header = header.split(None)
    pt_data = np.genfromtxt(pt_fl, skip_header=1).T

    pt_dict = {}
    for n in range(len(header)):
        pt_dict[header[n]] = pt_data[n, :]

    return pt_dict

def plot_fullspectrum(experiment_name, highres_wl, highres_flux,
                      lowres_wl, lowres_flux,
                      pt_dict):

    # fig, axis = plt.subplots(ncols=2, nrows=2)
    xtest = np.arange(0, 1, 0.01)
    ytest = np.sin(xtest)

    plt.figure()
    plt.title('title')

    ax1 = plt.subplot(221)
    ax1.plot(pt_dict['T'], pt_dict['P'])
    ax1.set_xlabel('Temperature [K]')
    ax1.set_ylabel('Pressure [bar]')
    ax1.set_yscale('log')
    ax1.invert_yaxis()
    ax1.set_title('Structure')

    keys = pt_dict.keys()
    keys.remove('P')
    keys.remove('T')

    ax2 = plt.subplot(222)
    N2 = np.ones_like(pt_dict['P'])
    for key in keys:
        ax2.plot(pt_dict[key], pt_dict['P'], label=key)
    ax2.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_ylabel('Pressure [bar]')
    ax2.set_xlabel('Mixing Ratio')
    ax2.set_title('Atmospheric Species')
    ax2.invert_yaxis()
    ax2.legend(bbox_to_anchor=(1, 0.9))


    ax3 = plt.subplot(212)


    ax3.set_xlabel(r'$\lambda$ [$\mu$m]')
    if args.spec_type == 'transit':
        #plt.ylabel(r'Transit Depth [(R$_{\mathrm{planet}}$ / R$_{\mathrm{star}}$)$^2$]')
        ax3.set_ylabel('Planet Radius [km]')
        spec = 'transit'
    elif args.spec_type =='direct':
        #plt.ylabel(r'Flux [W/m$^2$/$\mu$m]')
        ax3.set_ylabel(r'Reflectivity')
        spec = 'direct'

    ax3.plot(highres_wl, highres_flux, color='darkslategray',
            label='mixed, R=%i' % highres_R, linewidth=1)
    ax3.plot(lowres_wl, lowres_flux, color='paleturquoise',
            linestyle='dashed', label='mixed, R=%i' % lowres_R)

    ax3.legend()
    # #plt.show()


    # pos = ax3.get_position()
    # new_pos = [pos.x0, pos.y0, pos.x1, pos.y1]
    # ax3.set_position(new_pos)

    plt.tight_layout()



    plt.savefig('../plots/%s_%s_full.pdf' % (experiment_name, spec))

def extract_and_plot_full(experiment_name):
    if args.spec_type == 'transit':
        highres_fl = smart.readsmart.Trnst(trnst_file(experiment_name, 0.95, 0.55))


        # (r_planet/r_star)^2
#             highres_flux = highres_fl.tdepth

        # absolute radius
        highres_flux = highres_fl.absrad


    elif args.spec_type == 'direct':
        highres_fl = smart.readsmart.Rad(rad_file(experiment_name, 0.95, 0.55))

        highres_flux = highres_fl.pflux / highres_fl.sflux

    highres_wl = highres_fl.lam

    lowres_wl, lowres_flux = degrade_spectrum(highres_flux, highres_wl,
                                              0.95, R=lowres_R, span=0.55)



    pt_dict = extract_pt_dict(PT_DIR + experiment_name + '.pt')

    plot_fullspectrum(experiment_name, highres_wl, highres_flux,
                 lowres_wl, lowres_flux, pt_dict)

#####################################################

############# get experiment name ##############

if run_type == 'bands':
    experiment_name_lower = 'lower_o2inv%i_%s_%s_%s' % (int(100*o2_inv),
                                                                    n2scale,
                                                                    pt_shape,
                                                                    run_type)
    experiment_name_upper = 'upper_o2inv%i_%s_%s_%s' % (int(100*o2_inv),
                                                                    n2scale,
                                                                    pt_shape,
                                                                    run_type)
    experiment_name_mixed = 'mixed_o2inv%i_%s_%s_%s' % (int(100*o2_inv),
                                                                    n2scale,
                                                                    pt_shape,
                                                                    run_type)


    pressure, o2_high = get_pt_data(experiment_name_upper)
    pressure, o2_low = get_pt_data(experiment_name_lower)
    pressure, o2_mixed = get_pt_data(experiment_name_mixed)



    fig, axis = plt.subplots(nrows=4, ncols=5, sharex='col')

    axis[0, 0].set_title(r'O$_2$ Abundance') # , P$_{\mathrm{trop}}$ = %.01f' % args.trop_loc)
    axis[0, 1].set_title(r'$0.63 \mu m$')
    axis[0, 2].set_title(r'$0.68 \mu m$')
    axis[0, 3].set_title(r'$0.76 \mu m$')
    axis[0, 4].set_title(r'$1.27 \mu m$')

    axis[0, 0].loglog(o2_high, pressure, color='C0', linewidth=4)
    axis[1, 0].loglog(o2_low, pressure, color='C1', linewidth=4)
    axis[2, 0].loglog(o2_mixed, pressure, color='C2', linewidth=4)
    axis[3, 0].loglog(o2_high, pressure, label='upper', linewidth=4)
    axis[3, 0].loglog(o2_low, pressure, label='lower', linewidth=4)
    axis[3, 0].loglog(o2_mixed, pressure, label='mixed', linewidth=4)
    axis[3, 0].legend()

    for n in range(4):
        axis[n, 0].set_ylabel('P [bar]')
        axis[n, 0].set_xlabel(r'f(O$_2$)')
        axis[n, 0].invert_yaxis()

    bands= [0.63, 0.68, 0.76, 1.27]

    for n, band in enumerate(bands):

        highres_wl_high, highres_flux_high, lowres_wl_high, lowres_flux_high = extract(experiment_name_upper, band)
        highres_wl_low, highres_flux_low, lowres_wl_low, lowres_flux_low = extract(experiment_name_lower, band)
        highres_wl_mixed, highres_flux_mixed, lowres_wl_mixed, lowres_flux_mixed = extract(experiment_name_mixed, band)

        axis[0, n+1].plot(highres_wl_high, highres_flux_high, color='C0', label='upper, R=100,000')
        axis[0, n+1].plot(lowres_wl_high, lowres_flux_high, color='k', linewidth=2, linestyle='dashed', label='upper, R=150')

        axis[1, n+1].plot(highres_wl_low, highres_flux_low, color='C1', label='lower, R=100,000')
        axis[1, n+1].plot(lowres_wl_low, lowres_flux_low, color='k', linewidth=2, linestyle='dashed', label='lower, R=150')

        axis[2, n+1].plot(highres_wl_mixed, highres_flux_mixed, color='C2', label='mixed, R=100,000')
        axis[2, n+1].plot(lowres_wl_mixed, lowres_flux_mixed, color='k', linewidth=2, linestyle='dashed', label='mixed, R=150')


        axis[3, n+1].plot(highres_wl_mixed, highres_flux_mixed, color='C2', label='mixed, R=100,000')
        axis[3, n+1].plot(lowres_wl_mixed, lowres_flux_mixed, color='C2', linewidth=3, linestyle='dashed', label='mixed, R=150', alpha=0.7)
        axis[3, n+1].plot(highres_wl_low, highres_flux_low, color='C1', label='lower, R=100,000')
        axis[3, n+1].plot(lowres_wl_low, lowres_flux_low, color='C1', linewidth=3, linestyle='dashed', label='lower, R=150', alpha=0.7)
        axis[3, n+1].plot(highres_wl_high, highres_flux_high, color='C0', label='upper, R=100,000')
        axis[3, n+1].plot(lowres_wl_high, lowres_flux_high, color='C0', linewidth=3, linestyle='dashed', label='upper, R=150', alpha=0.7)

        axis[0, 4].legend(bbox_to_anchor=(1, 0.9))
        axis[1, 4].legend(bbox_to_anchor=(1, 0.9))
        axis[2, 4].legend(bbox_to_anchor=(1, 0.9))
        axis[3, 4].legend(bbox_to_anchor=(1, 0.9))



    for n in range(4):
        for m in range(5):
            axis[n, m].tick_params(axis='both', reset=True)



    if args.spec_type == 'direct':
        spec = 'direct'
        for n in range(4):
            for m in range(4):
                axis[n, m+1].set_ylabel('Reflectivity')
                axis[n, m+1].set_xlabel(r'$\lambda$ [$\mu m$]')

    elif args.spec_type == 'transit':
        spec = 'transit'
        for n in range(4):
            for m in range(4):
                axis[n, m+1].set_ylabel('Planet Radius [km]')
                axis[n, m+1].set_xlabel(r'$\lambda$ [$\mu m$]')
    xtxt = -0.01

    #fig.text(0.5, 1.01, args.title, fontsize=40, fontweight='bold', va='center', ha='center')
    fig.text(xtxt, 0.84, 'Upper', fontsize=25, fontweight='bold', rotation=90, va='center', ha='center')
    fig.text(xtxt, 0.62, 'Lower', fontsize=25, fontweight='bold', rotation=90, va='center', ha='center')
    fig.text(xtxt, 0.40, 'Mixed', fontsize=25, fontweight='bold', rotation=90, va='center', ha='center')
    fig.text(xtxt, 0.16, 'All', fontsize=25, fontweight='bold', rotation=90, va='center', ha='center')

    fig.tight_layout()

    plt.savefig('../plots/%s_%s_report.pdf' % (experiment, spec), bbox_inches='tight')


elif run_type == 'full':
    extract_and_plot_full(experiment)
