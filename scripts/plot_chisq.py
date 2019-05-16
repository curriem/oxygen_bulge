import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import smart

spectrum_type = 'transit'
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

min_ind = np.where(full_chisq == np.min(full_chisq))
optimal_o2_val = full_o2_vals[min_ind]

plt.figure()
plt.subplot(521)

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

plt.subplot(522)
for n in range(len(full_o2_vals)):
    plt.scatter(full_o2_vals[n], full_chisq[n], color='k')

plt.scatter(optimal_o2_val, np.min(full_chisq), s=200, marker='*', color='C0', label='optimum mixed o2')

plt.xlabel(r'f(O$_2$) [mol/mol]')
plt.ylabel(r'$\chi^2$')
plt.xlim(0, .2)
plt.legend()


plt.subplot(512)
plt.plot(wl, observation_flux, label='upper o2', color='C2')

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

plt.plot(wl, optimal_model_flux, label='mixed o2', color='C0')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(ylab)
plt.legend()

plt.subplot(513)
plt.plot(wl, observation_flux - optimal_model_flux, color='k')
plt.xlabel(r'$\lambda$ [$\mu m$]')
plt.ylabel(r'Upper O$_2$ $-$ Mixed O$_2$')

wlrange = 0.04
bands = [0.63, 0.68, 0.76, 1.27]
for n, band in enumerate(bands):
    plt.subplot(5, 4, n+13)
    wlmax = band + wlrange
    wlmin = band - wlrange
    band_inds = (wl > wlmin) & (wl < wlmax)
    plt.plot(wl[band_inds], observation_flux[band_inds], color='C2')
    plt.plot(wl[band_inds], optimal_model_flux[band_inds], color='C0')
    plt.title(r'%s $\mu m$ O$_2$ band' % str(band))
    plt.xlabel(r'$\lambda$ [$\mu m$]')
    plt.ylabel(ylab)


for n, band in enumerate(bands):
    plt.subplot(5, 4, n+17)
    wlmax = band + wlrange
    wlmin = band - wlrange
    band_inds = (wl > wlmin) & (wl < wlmax)
    plt.plot(wl[band_inds], observation_flux[band_inds] - optimal_model_flux[band_inds], color='k')
    plt.title(r'Upper O$_2$ $-$ Mixed O$_2$')
    plt.xlabel(r'$\lambda$ [$\mu m$]')

plt.tight_layout()
plt.savefig('%s_degeneracy_report.pdf' % spectrum_type)
