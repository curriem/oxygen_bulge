import numpy as np
import smart
import project_tools
import glob
import sys

spectrum_type = sys.argv[1]

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

newrange = (wl > 0.4)

wl = wl[newrange]
observation_flux = observation_flux[newrange]



chisqs = []
o2_vals = []
for model_fl in models:
    if spectrum_type == 'direct':
        model = smart.readsmart.Rad(model_fl)
        model_flux = model.pflux / model.sflux
    elif spectrum_type == 'transit':
        model = smart.readsmart.Trnst(model_fl)
        model_flux = model.absrad
    model_flux = model_flux[new_range]
    chisq = calculate_chisq(observation_flux, model_flux)
    chisqs.append(chisq)

    o2_val = '0.' + model_fl.split('_')[o2val_ind]
    o2_val = float(o2_val)
    o2_vals.append(o2_val)
    print 'O2:', o2_val
    print 'Chi^2:', chisq
    print '\n'

np.savetxt('%s_trappist1_t1e_degeneracy_metadata.txt' % spectrum_type,np.array([o2_vals, chisqs]))
