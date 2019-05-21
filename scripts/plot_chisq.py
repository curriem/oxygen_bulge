import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import smart
import sys

spectrum_type = sys.argv[1]
####### fix plotting default values #######

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

figure = {'figsize' : (14, 16)}

lines = {'linewidth' : 1}

matplotlib.rc('font', **font)
matplotlib.rc('figure', **figure)
matplotlib.rc('lines', **lines)

#############################################

PT_DIR = '../data/pt_fls/'
SMART_OUTPUTS = '../data/smart_outputs/'

full_fl = np.loadtxt('%s_trappist1_t1e_degeneracy_metadata.txt' % spectrum_type)
bands_fl = np.loadtxt('%s_trappist1_t1e_degeneracy_metadata_bands.txt' % spectrum_type).T

full_o2_vals, full_chisq = full_fl[0], full_fl[1]

bands_o2_vals, bands_band, bands_chisq = bands_fl[0], bands_fl[1], bands_fl[2]

observation_pt = np.genfromtxt(PT_DIR + 'upper_oxygen_0_2.pt', skip_header=1)
observation_p = observation_pt[:, 0]
observation_o2 = observation_pt[:, 8]

if spectrum_type == 'direct':
    observation = smart.readsmart.Rad(SMART_OUTPUTS + 'trappist1_t1e_upper_oxygen_0_2_6666_250000cm_toa.rad')
    observation_flux = observation.pflux / observation.sflux
elif spectrum_type == 'transit':
    observation = smart.readsmart.Trnst(SMART_OUTPUTS + 'trappist1_t1e_upper_oxygen_0_2_6666_250000cm.trnst')
    observation_flux = observation.absrad
wl = observation.lam

newrange = (wl > 0.4)

wl = wl[newrange]
observation_flux = observation_flux[newrange]

min_ind = np.where(full_chisq == np.min(full_chisq))
optimal_o2_val = full_o2_vals[min_ind]

plt.figure()
plt.subplot(621)

for o2_val in full_o2_vals:
    plt.axvline(o2_val, linewidth=0.5, color='k')
plt.plot(observation_o2, observation_p, color='C2', linewidth=5, label='upper o2')
plt.axvline(optimal_o2_val, color='C0', linewidth=5, label='mixed o2')
plt.gca().invert_yaxis()
plt.yscale('log')
plt.ylim(1, 1e-7)
plt.xlabel(r'f(O$_2$) [mol/mol]')
plt.ylabel(r'p [bar]')
plt.legend()

plt.subplot(622)
for n in range(len(full_o2_vals)):
    plt.scatter(full_o2_vals[n], full_chisq[n], color='k')

plt.scatter(optimal_o2_val, np.min(full_chisq), s=200, marker='*', color='C0', label='optimum mixed o2')

plt.xlabel(r'f(O$_2$) [mol/mol]')
plt.ylabel(r'$\chi^2$')
plt.legend()


plt.subplot(612)
plt.plot(wl, observation_flux, label='upper o2: 0.2', color='k')

if spectrum_type == 'direct':
    optimal_model_fl = SMART_OUTPUTS + 'trappist1_t1e_mixed_oxygen_0_%s_6666_250000cm_toa.rad' % str(optimal_o2_val[0]).split('.')[-1]
    optimal_model = smart.readsmart.Rad(optimal_model_fl)
    optimal_model_flux = optimal_model.pflux / optimal_model.sflux
    ylab = 'Reflectivity'

elif spectrum_type == 'transit':
    optimal_model_fl = SMART_OUTPUTS + 'trappist1_t1e_mixed_oxygen_0_%s_6666_250000cm.trnst' % str(optimal_o2_val[0]).split('.')[-1]
    optimal_model = smart.readsmart.Trnst(optimal_model_fl)
    optimal_model_flux = optimal_model.absrad
    ylab = 'Absolute Radius [km]'
wl = observation.lam

overall_wl = wl[newrange]
optimal_model_flux = optimal_model_flux[newrange]

plt.plot(overall_wl, optimal_model_flux, label='mixed o2: %s' % str(optimal_o2_val[0]), color='C0', alpha=0.5)
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(ylab)
plt.legend()

plt.subplot(613)
plt.plot(overall_wl, observation_flux - optimal_model_flux, color='k')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'Upper O$_2$ $-$ Mixed O$_2$')

wlrange = 0.04
bands = [0.63, 0.68, 0.76, 1.27]

for n, band in enumerate(bands):

    band_inds = np.where(bands_band == band)
    plt.subplot(6, 4, n+13)
    plt.scatter(bands_o2_vals[band_inds], bands_chisq[band_inds], color='k')
    band_optimal_o2_val = bands_o2_vals[band_inds][np.argmin(bands_chisq[band_inds])]
    plt.scatter(band_optimal_o2_val, np.min(bands_chisq[band_inds]), s=200, marker='*', color='C0', label='optimum mixed o2')
    plt.xlabel(r'f(O$_2$) [mol/mol]')
    plt.ylabel(r'$\chi^2$')
    plt.legend()

    if spectrum_type == 'direct':
        band_optimal_model_fl = SMART_OUTPUTS + 'trappist1_t1e_mixed_oxygen_0_%s_6666_250000cm_toa.rad' % str(band_optimal_o2_val).split('.')[-1]
        print band_optimal_model_fl
        band_optimal_model = smart.readsmart.Rad(band_optimal_model_fl)
        band_optimal_model_flux = band_optimal_model.pflux / band_optimal_model.sflux

    elif spectrum_type == 'transit':
        band_optimal_model_fl = SMART_OUTPUTS + 'trappist1_t1e_mixed_oxygen_0_%s_6666_250000cm.trnst' % str(band_optimal_o2_val).split('.')[-1]
        print band_optimal_model_fl
        band_optimal_model = smart.readsmart.Trnst(band_optimal_model_fl)
        band_optimal_model_flux = band_optimal_model.absrad

    band_optimal_model_flux = band_optimal_model_flux[newrange]
    plt.subplot(6, 4, n+17)


    wlmax = band + wlrange
    wlmin = band - wlrange
    band_inds = (overall_wl > wlmin) & (overall_wl < wlmax)

    plt.plot(overall_wl[band_inds], observation_flux[band_inds], color='k')
    plt.plot(overall_wl[band_inds], band_optimal_model_flux[band_inds], color='C0', label='o2: %s' % str(band_optimal_o2_val), alpha=0.5)
    plt.title(r'%s $\mu m$ O$_2$ band' % str(band))
    plt.xlabel(r'$\lambda$ [$\mu m$]')
    plt.ylabel(ylab)
    plt.legend()


    plt.subplot(6, 4, n+21)

    plt.plot(overall_wl[band_inds], observation_flux[band_inds] - band_optimal_model_flux[band_inds], color='k')
    plt.title(r'Upper O$_2$ $-$ Mixed O$_2$')
    plt.xlabel(r'$\lambda$ [$\mu m$]')

plt.tight_layout()
plt.savefig('%s_degeneracy_report.pdf' % spectrum_type)
