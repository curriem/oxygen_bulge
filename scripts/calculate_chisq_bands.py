import numpy as np
import smart
import project_tools
import glob
spectrum_type = 'transit'
SMART_OUTPUTS = '../data/smart_outputs/'
PT_DIR = '../data/pt_fls/'

def calculate_chisq(observation, model):
    return np.sum((observation - model)**2)

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


metadata = []
bands = [0.63, 0.68, 0.76, 1.27]
wlrange = 0.04
for model_fl in models:
    if spectrum_type == 'direct':
        model = smart.readsmart.Rad(model_fl)
        model_flux = model.pflux / model.sflux
    elif spectrum_type == 'transit':
        model = smart.readsmart.Trnst(model_fl)
        model_flux = model.absrad
    o2_val = '0.' + model_fl.split('_')[o2val_ind]
    o2_val = float(o2_val)
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
        metadata.append([o2_val, band, chisq])

np.savetxt('%s_trappist1_t1e_degeneracy_metadata_bands.txt' % spectrum_type,np.array(metadata))
