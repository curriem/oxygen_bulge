import os
import stat
import smart
import numpy as np
try:
    from scipy.interpolate import interp1d
except:
    pass


PT_DIR = '../data/pt_fls/'

def run_trappist(ptfile, band, wlrange, place='output', R=100000., rm_when_done=False,
              o2o2=True, transit=True, NCPU=1, addn2=True):

        '''
        Parameters
        -----------
        place : str
            Directory to create to run this test in
        rm_when_done : bool
            Set to remove ``place'' when the simulation finishes
        R : float
            resolving power
        '''
        try:
            os.mkdir(place)
        except OSError:
            pass

        # load pt file (relative import frim photochem directory)
        atm = smart.interface.Atmosphere.from_pt(ptfile, atm_dir=place,
                                                 addn2=addn2)

        # Set atmosphere parameters
        atm.set_dataframe_from_profiles()
        atm.write_atm()

        # Prep simulation
        sim = smart.interface.Smart(tag='case', atmosphere=atm)
        sim.smartin.source = 3 # solar and thermal sources
        sim.smartin.out_format = 1 # no transit calculations

        sim.smartin.spec = '../../specs/TRAPPIST-1.dat'
        sim.smartin.spec_skip = 5
        sim.smartin.spec_unit = 2 # Set unit to W/m^2/um

        sim.smartin.r_AU = 0.02916

        sim.smartin.radius = 0.915 * 6371.
        sim.smartin.grav = 0.930 * 9.8

        sim.lblin.radius = 0.915 * 6371.
        sim.lblin.grav = 0.930 * 9.8

        sim.set_run_in_place(place=place)

        wlmin = band - wlrange
        wlmax = band + wlrange

        # change the wavelength range
        sim.smartin.minwn = 1e4 / wlmax
        sim.smartin.maxwn = 1e4 / wlmin
        sim.lblin.minwn = 1e4 / wlmax
        sim.lblin.maxwn = 1e4 / wlmin

        # watch out for indexing: check if name equals the gas I want
        if not o2o2:
            sim.atmosphere.gases[6].cia_file = None

        if transit:
            # change to transit:
            sim.smartin.out_format = 9
            sim.smartin.itrnst = 2 # ray tracing
            sim.smartin.radstar = 0.121 # times sun radius
            sim.smartin.irefract = 1 # turn on refraction


        # modify resolution of spectrum
        delt_nu = 1e4 / ((wlmax - wlmin)/2 * R)
        sim.smartin.FWHM = delt_nu
        sim.smartin.sample_res = delt_nu

        # set to run on this machine
        sim.set_executables_automatically()

        # generate LBLABC scripts
        sim.gen_lblscripts()

        # run LBLABC on all but two CPUs on machine
        sim.run_lblabc(ncpu=NCPU)

        #sim.tag = sim.tag + '_' + tag + '_' + 'R' + str(R)

        # run SMART
        sim.write_smart(write_file=True)
        sim.run_smart()

        # remove temp directory
        if rm_when_done:
            os.system('rm -r %s' % place)

        return


def run_sun(ptfile, band, wlrange, place='output', R=100000., rm_when_done=False,
              o2o2=True, transit=True, NCPU=1, addn2=True):

        '''
        Parameters
        -----------
        place : str
            Directory to create to run this test in
        rm_when_done : bool
            Set to remove ``place'' when the simulation finishes
        R : float
            resolving power
        '''
        try:
            os.mkdir(place)
        except OSError:
            pass

        # load pt file (relative import frim photochem directory)
        atm = smart.interface.Atmosphere.from_pt(ptfile, atm_dir=place,
                                                 addn2=addn2)

        # Set atmosphere parameters
        atm.set_dataframe_from_profiles()
        atm.write_atm()

        # Prep simulation
        sim = smart.interface.Smart(tag='case', atmosphere=atm)
        sim.smartin.source = 3 # solar and thermal sources
        sim.smartin.out_format = 1 # no transit calculations


        sim.set_run_in_place(place=place)

        wlmin = band - wlrange
        wlmax = band + wlrange

        # change the wavelength range
        sim.smartin.minwn = 1e4 / wlmax
        sim.smartin.maxwn = 1e4 / wlmin
        sim.lblin.minwn = 1e4 / wlmax
        sim.lblin.maxwn = 1e4 / wlmin

        # watch out for indexing: check if name equals the gas I want
        if not o2o2:
            sim.atmosphere.gases[6].cia_file = None

        if transit:
            # change to transit:
            sim.smartin.out_format = 9
            sim.smartin.itrnst = 2 # ray tracing
            sim.smartin.radstar = 1 # times sun radius
            sim.smartin.irefract = 1 # turn on refraction


        # modify resolution of spectrum
        delt_nu = 1e4 / ((wlmax - wlmin)/2 * R)
        sim.smartin.FWHM = delt_nu
        sim.smartin.sample_res = delt_nu

        # set to run on this machine
        sim.set_executables_automatically()

        # generate LBLABC scripts
        sim.gen_lblscripts()

        # run LBLABC on all but two CPUs on machine
        sim.run_lblabc(ncpu=NCPU)

        #sim.tag = sim.tag + '_' + tag + '_' + 'R' + str(R)

        # run SMART
        sim.write_smart(write_file=True)
        sim.run_smart()

        # remove temp directory
        if rm_when_done:
            os.system('rm -r %s' % place)

        return


def make_pt_fl(max_o2_abundance, o2_location, trop_loc, h2o, pt_shape, experiment_name,
                N2_scalar=1):
    # Standard Earth O_2 mixing ratio:
    earth_standard_ptfile = PT_DIR + '/earth_standard_icrccm_vmix.pt'
    header = ["Press",
              "Temp",
              "H2O",
              "CO2",
              "O3",
              "N2O",
              "CO",
              "CH4",
              "O2",
              "N2"]

    header = "        ".join(header)
    earth_pt_data = np.genfromtxt(earth_standard_ptfile, skip_header=1)
    earth_pt_data_inv = earth_pt_data.T
    earth_pressure = earth_pt_data_inv[0, :]
    earth_o2_mixing_ratio = earth_pt_data_inv[8, :]

    #new_pt_fl_name = 'o2%s%s_water%i.pt' % (o2_location, str(int(max_o2_abundance*100)).zfill(2), int(h2o))
    new_pt_fl_name = experiment_name + '.pt'
    if pt_shape == 'wedge':
        pressure = [0, np.log10(trop_loc), -7]
        #pressure = [0, -1, -7]
        if o2_location == 'lower':
            o2_mix = [np.log10(max_o2_abundance), -5, -14]
            fill_val = np.log10(max_o2_abundance)
        elif o2_location == 'upper':
            o2_mix = [-14, -5, np.log10(max_o2_abundance)]
            fill_val = -14
        else:
            # mixed o2 abundance
            o2_mix = [np.log10(max_o2_abundance)]*3
            fill_val = np.log10(max_o2_abundance)

    elif pt_shape == 'step':
        pressure = [0, np.log10(trop_loc), np.log10(trop_loc - 0.00001), -7]
        if o2_location == 'lower':
            o2_mix = [np.log10(max_o2_abundance), np.log10(max_o2_abundance), -14, -14]
            fill_val = np.log10(max_o2_abundance)
        elif o2_location == 'upper':
            o2_mix = [-14, -14, np.log10(max_o2_abundance), np.log10(max_o2_abundance)]
            fill_val = -14
        else:
            # mixed o2 abundance
            o2_mix = [np.log10(max_o2_abundance)]*4
            fill_val = np.log10(max_o2_abundance)


    func = interp1d(pressure, o2_mix, bounds_error=False,
                    fill_value=fill_val)
    new_o2_mix = func(np.log10(earth_pressure))
    new_o2_mix = 10**new_o2_mix

    new_pt_inv = np.copy(earth_pt_data_inv)
    new_pt_inv[8, :] = new_o2_mix

    if not h2o:
        new_pt_inv[2, :] /= 1000


    sum_mixing_ratios = np.sum(new_pt_inv[2:, :], axis=0)

    sum_pressures = earth_pressure * sum_mixing_ratios


    N2_earth_mixing_ratio = np.ones_like(sum_mixing_ratios) - sum_mixing_ratios

    N2_earth_pressure = N2_earth_mixing_ratio * earth_pressure

    N2_new_pressure = N2_scalar * N2_earth_pressure

    sum_new_pressures = sum_pressures + N2_new_pressure

    N2_new_mixing_ratio = N2_new_pressure / sum_new_pressures

    #print N2_pressure

    new_pt_inv = np.vstack((new_pt_inv, N2_new_mixing_ratio))

    new_pt_inv[0, :] = sum_new_pressures

    new_pt_data = new_pt_inv.T
    np.savetxt(PT_DIR+new_pt_fl_name, new_pt_data, delimiter='   ', header=header, comments='')

    return new_pt_fl_name

def plot_pt(pt_fl):
    pt = np.genfromtxt(PT_DIR + pt_fl, skip_header=1)
    pt_data = pt.T
    pressure = pt_data[0, :]
    temperature = pt_data[1, :]
    o2 = pt_data[8, :]
    plt.figure(figsize=(12, 10))
    plt.loglog(o2, pressure)
    plt.gca().invert_yaxis()
    plt.xlabel(r'O$_2$ mixing ratio')
    plt.ylabel('Pressure [bar]')
    plt.show()


def write_slurm_script_python(runfiles, name="smart_run", subname="smart_submit.csh",
                       workdir = "",
                       nodes = 1, mem = "500G", walltime = "0", ntasks = 28,
                       account = "vsm", submit = False, rm_after_submit = False,
                       preamble = ['module load icc_18',
                                   'module load parallel-20170722']):
    """
    Write a hyak SLURM bash script for simple python scripts

    Parameters
    ----------
    runfiles : array of strs
        Array of all python scripts to run
    name : str
        Name of the job
    subname : str
        Name of the bash submission script to generate
    workdir : str
        Working directory
    nodes : int
        Number of nodes
    mem : str
        Requested memory
    walltime : str
        Walltime to allow simulation to run for ("0" is indefinite)
    ntasks : int
        Number of processors to allow usage of
    account : str
        Name of account and partition on hyak
    submit : bool
        Set to submit the job
    rm_after_submit : bool
        Set to remove script immediately after submission
    preamble : list
        Bash statements to call before running python script
        (e.g. ['module load icc_18', 'source activate my_root'])
    """

    # Get absolute path of workdir
    abs_place = os.path.abspath(workdir)

    newfile = os.path.join(abs_place, subname)

    f = open(newfile, 'w')

    f.write('#!/bin/bash\n')
    f.write('\n')
    f.write('#SBATCH --job-name=%s\n' %name)
    f.write('\n')
    f.write('#SBATCH --account=%s\n' %account)
    f.write('#SBATCH --partition=%s\n' %account)
    f.write('\n')
    f.write('#SBATCH --nodes=%i\n' %nodes)
    ###SBATCH --ntasks-per-node=28 ### not necessarily required, so I'm not using it in this one
    f.write('#SBATCH --time=%s\n' %walltime)  ### this basically allows it to run for 365 days. Note it seems to otherwise only accept hours:mins:sec
    f.write('#SBATCH --mem=%s\n' %mem) ### apparently you only get the memory you ask for
    f.write('\n')
    f.write('#SBATCH --ntasks=%i\n' %ntasks)
    f.write('#SBATCH --exclusive\n')
    f.write('\n')
    f.write('#SBATCH --workdir=%s\n' %abs_place)
    f.write('\n')
    f.write('ulimit -s unlimited\n')

    # Loop over preamble statements
    for i in range(len(preamble)):
        f.write('%s\n' %preamble[i])

    f.write('\n')

    for runfile in runfiles:
        f.write('python %s\n' %runfile)

    f.close()

    # Set permissions to allow execute
    st = os.stat(newfile)
    os.chmod(newfile, st.st_mode | stat.S_IEXEC)

    # Submit the run to the queue?
    if submit:

        subprocess.check_call(["sbatch", "-p", "vsm", "-A", "vsm", newfile])

        # Delete the bash script after submitting
        if rm_after_submit:

            os.system('rm %s' %newfile)

    return
