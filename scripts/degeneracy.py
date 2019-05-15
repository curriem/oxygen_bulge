import smart
import numpy as np
import argparse
import matplotlib.pyplot as plt
import project_tools
import glob
import matplotlib

'''
this script attempts to vary O2 abundance in a well mixed atmosphere to match
the O2 spectrum in an oxygen bulge scenario for a give band +/- 0.04 microns
'''
def extract_pt_dict(pt_fl):

    header = "P        T       H2O             CO2        O3            N2O          CO           CH4         O2     N2"
    header = header.split(None)
    pt_data = np.genfromtxt(pt_fl, skip_header=1).T

    pt_dict = {}
    for n in range(len(header)):
        pt_dict[header[n]] = pt_data[n, :]

    return pt_dict


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}

figure = {'figsize' : (10, 4)}

lines = {'linewidth' : 2}

matplotlib.rc('font', **font)
matplotlib.rc('figure', **figure)
matplotlib.rc('lines', **lines)


parser = argparse.ArgumentParser()
parser.add_argument('-system', dest='system', type=str)
args = parser.parse_args()

upper_bulge_experiment = '../data/smart_outputs/%s/upper_o2inv21_n2scale1_step_bands_output/case_12500_13888cm_toa.rad' % args.system
o2_abundance_experiments = glob.glob('../data/smart_outputs/degeneracy/%s/*/case_12500_13888cm_toa.rad' % args.system)
o2_abundance_experiments = sorted(o2_abundance_experiments)

upper_bulge_pt = '../data/pt_fls/upper_o2inv21_n2scale1_step_bands.pt'
o2_abundance_pts = glob.glob('../data/pt_fls/*single*')
o2_abundance_pts = sorted(o2_abundance_pts)

plt.figure()
plt.subplot(121)
plt.ylabel(r'Pressure [bar]')
plt.xlabel(r'f(O$_2$)')
pt_dict = extract_pt_dict(upper_bulge_pt)
plt.plot(pt_dict['O2'], pt_dict['P'], color='k', linewidth=5, alpha=0.6)

plt.subplot(122)
plt.xlabel(r'$\lambda$ ($\mu m$)')
plt.ylabel(r'Reflectivity')
upper_bulge_open = smart.readsmart.Rad(upper_bulge_experiment)
upper_bulge_flux = upper_bulge_open.pflux / upper_bulge_open.sflux
upper_bulge_wl = upper_bulge_open.lam

plt.plot(upper_bulge_wl, upper_bulge_flux, color='k', linewidth=5, alpha=0.6)

for n, exp in enumerate(o2_abundance_experiments[:5]):
    plt.subplot(121)
    pt_dict = extract_pt_dict(o2_abundance_pts[n])
    plt.plot(pt_dict['O2'], pt_dict['P'], alpha=0.3)

    plt.subplot(122)
    exp_open = smart.readsmart.Rad(exp)
    exp_flux = exp_open.pflux / exp_open.sflux
    exp_wl = exp_open.lam
    plt.plot(exp_wl, exp_flux, alpha=0.3)

plt.tight_layout()
plt.savefig('degeneracy.pdf')
