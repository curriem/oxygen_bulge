import numpy as np
import smart
import project_tools
import glob

SMART_OUTPUTS = '../data/smart_outputs/'
PT_DIR = '../data/pt_fls/'

def calculate_chisq(observation, model):
    return np.sum((observation - model)**2)

observation = smart.readsmart.Rad(SMART_OUTPUTS + 'trappist1_t1e_upper_oxygen_0_2_6666_250000cm_toa.rad')
observation_flux = observation.pflux / observation.sflux
wl = observation.lam

models = glob.glob(SMART_OUTPUTS+ 'trappist1_t1e_mixed_oxygen_0_*_6666_250000cm_toa.rad')


metadata = []
bands = [0.63, 0.68, 0.76, 1.27]
wlrange = 0.04
for model in models:
    model = smart.readsmart.Rad(model_fl)
    model_flux = model.pflux / model.sflux
    o2_val = '0.' + model_fl.split('_')[-4]
    o2_val = float(o2_val)
    for band in bands:

        wlmax = band + wlrange
        wlmin = band - wlrange

        band_inds = (wl > wlmin) & (wl < wlmax)
        band_wl = wl[band_inds]
        band_obs_flux = observation_flux[band_inds]

        band_model_flux = model_flux[band_inds]

        chisq = calculate_chisq(band_obs_flux, band_model_flux)

        metadata.append([o2_val, band, chisq])

np.savetxt('trappist1_t1e_degeneracy_metadata_bands.txt',[o2_vals, chisqs])
