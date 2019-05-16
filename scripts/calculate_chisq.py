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



observation_pt = np.genfromtxt(PT_DIR + 'upper_oxygen_0_2.pt', skip_header=1)
observation_p = observation_pt[:, 0]
observation_o2 = observation_pt[:, 8]


models = glob.glob(SMART_OUTPUTS+ 'trappist1_t1e_mixed_oxygen_0_*_6666_250000cm_toa.rad')
chisqs = []
o2_vals = []
for model_fl in models:
    model = smart.readsmart.Rad(model_fl)
    model_flux = model.pflux / model.sflux
    assert np.array_equal(model.lam, wl)
    chisq = calculate_chisq(observation_flux, model_flux)
    chisqs.append(chisq)

    o2_val = '0.' + model_fl.split('_')[-4]
    o2_val = float(o2_val)
    o2_vals.append(o2_val)
    print o2_val
    print chisq

np.savetxt('trappist1_t1e_degeneracy_metadata.txt',[o2_vals, chisqs])
