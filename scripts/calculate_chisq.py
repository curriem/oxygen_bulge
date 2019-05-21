import numpy as np
import smart
import project_tools
import glob
import sys
import coronagraph as cg

R = 100000
if R != 100000:
    coro = True
else:
    coro = False


spectrum_type = sys.argv[1]

SMART_OUTPUTS = '../data/smart_outputs/'
PT_DIR = '../data/pt_fls/'

def calculate_chisq(observation, model):
    return np.sum((observation - model)**2)

def reduce_R(old_wl, old_flux):
    new_wl, dlam = cg.noise_routines.construct_lam(np.min(old_wl), np.max(old_wl), Res=R)
    new_flux = cg.downbin_spec(old_flux, new_wl, old_wl, dlam=dlam)
    return new_wl, new_flux

if spectrum_type == 'direct':
    observation = smart.readsmart.Rad(SMART_OUTPUTS + 'trappist1_t1e_upper_oxygen_0_2_6666_250000cm_toa.rad')
    observation_flux = observation.pflux / observation.sflux
    models = glob.glob(SMART_OUTPUTS+ 'trappist1_t1e_mixed_oxygen_0_*_6666_250000cm_toa.rad')
    o2val_ind = -4

elif spectrum_type == 'transit':
    observation = smart.readsmart.Trnst(SMART_OUTPUTS + 'trappist1_t1e_upper_oxygen_0_2_6666_250000cm.trnst')
    observation_flux = observation.absrad
    models = glob.glob(SMART_OUTPUTS+ 'trappist1_t1e_mixed_oxygen_0_*_6666_250000cm.trnst')
    o2val_ind = -3



wl = observation.lam

newrange = (wl > 0.4)

wl = wl[newrange]
observation_flux = observation_flux[newrange]

if coro:
    wl, observaition_flux = reduce_R(wl, observation_flux)

metadata_full = []
metadata_bands = []
wlrange = 0.04
bands = [0.63, 0.68, 0.76, 1.27]
for model_fl in models:
    if spectrum_type == 'direct':
        model = smart.readsmart.Rad(model_fl)
        model_flux = model.pflux / model.sflux
    elif spectrum_type == 'transit':
        model = smart.readsmart.Trnst(model_fl)
        model_flux = model.absrad
    model_wl = model.lam[newrange]
    model_flux = model_flux[newrange]
    if coro:
        model_wl, model_flux = reduce_R(model_wl, model_flux)

    chisq = calculate_chisq(observation_flux, model_flux)
    o2_val = '0.' + model_fl.split('_')[o2val_ind]
    o2_val = float(o2_val)
    metadata_full.append([o2_val, chisq])
    print 'O2:', o2_val
    print 'Chi^2:', chisq
    print '\n'

    for band in bands:

        wlmax = band + wlrange
        wlmin = band - wlrange

        band_inds = (wl > wlmin) & (wl < wlmax)
        band_wl = wl[band_inds]
        band_obs_flux = observation_flux[band_inds]

        band_model_flux = model_flux[band_inds]

        chisq = calculate_chisq(band_obs_flux, band_model_flux)

        print 'O2:', o2_val
        print 'Band:', band
        print 'Chi^2:', chisq
        print '\n'
        metadata_bands.append([o2_val, band, chisq])

np.savetxt('%s_trappist1_t1e_degeneracy_metadata_R%i.txt' % (spectrum_type, int(R)),np.array(metadata_full))
np.savetxt('%s_trappist1_t1e_degeneracy_metadata_bands_R%i.txt' % (spectrum_type, int(R)) ,np.array(metadata_bands))
