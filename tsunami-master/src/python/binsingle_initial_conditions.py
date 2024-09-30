"""
binsingle_initial_conditions.py

This script generates initial conditions for a set of binary-single hyperbolic encounters.
It requires pandas and scipy installed, including hdf5 libraries.
The output file is a pandas hdf5 file, named <setname>.h5, which contains the initial conditions in tabular data.
Each row of this file is a single initial condition for an individual binary-single encounter.

This script also creates a bash file <setname>.sh, which creates a list of individual runs and
runs it in parallel using gnu parallel. Each individual run is executed via the script binsingle_run.py

The function initial_conditions_isol is the function that defines and samples the initial condition distributions
of the three-body encounters, like impact parameter, binary properties, velocity at infinity.

You can create your own function to define the initial conditions similar to initial_conditions_isol.
Just create your function and pass it as an arguments to generate_initial_condition_set,
together with the function arguments in a dict and passed to ic_args, as in:
ic_func = <your_function_name>
ic_args = <{"dictionary" : with_your_args}>

Example usage:
    * python binsingle_initial_conditions.py testset -N 100
    * bash testset.sh

Creates a set named testset, with 100 binary-single encounters, using the default initial condition parameters
Check the varius input parameters in the argument options with -h
"""


import tsunami
import numpy as np
import os
import shutil
import pandas as pd
from scipy.stats import loguniform
from scipy.stats import maxwell
from binsingle_run import output_suff, coll_folder_suff
from binsingle_evolve_system import run_options
import argparse

KU = tsunami.KeplerUtils()


def initial_conditions_isol(Mmin=10, Mmax=50, bmax_ainn=2, ainn_min=0.1, ainn_max=1, sigma_vinf=1.0, Rschw_mult=10):
    m1 = loguniform.rvs(Mmin, Mmax)
    m2 = loguniform.rvs(Mmin, Mmax)
    m3 = loguniform.rvs(Mmin, Mmax)

    # Set up collision radii, Schwarchild radius for BHs, tidal radius for star (G=1)
    R1 = Rschw_mult * 2 * m1 / KU.speed_of_light ** 2
    R2 = Rschw_mult * 2 * m2 / KU.speed_of_light ** 2
    R3 = Rschw_mult * 2 * m3 / KU.speed_of_light ** 2

    a_inn = np.random.uniform(ainn_min, ainn_max)
    b = a_inn * bmax_ainn * np.random.uniform(0.0, 1.0)**0.5
    e_inn = 0.0
    ome_inn = 0.0

    v_to_kms = 2 * np.pi * KU.au2km / KU.yr2sec
    sigma_vinf_nb = sigma_vinf/v_to_kms
 #   while True:
    vinf = maxwell.rvs(scale=sigma_vinf_nb)
#        if vinf < 2*sigma_vinf_nb: break

    M_inn = np.random.uniform(0, 2 * np.pi)
    nu_inn = KU.M_to_nu(M_inn, e_inn)

    i_out = np.random.uniform(-1, 1)
    i_out = np.arccos(i_out)
    Ome_out = np.random.uniform(0, 2 * np.pi)
    ome_out = np.random.uniform(0, 2 * np.pi)

    data = {"m1": m1, "m2": m2, "m3": m3, "R1": R1, "R2": R2, "R3": R3,
            "a_inn": a_inn, "e_inn": e_inn, "nu_inn": nu_inn, "ome_inn": ome_inn,
            "vinf": vinf, "b": b,  "i_out": i_out, "ome_out": ome_out, "Ome_out": Ome_out}

    return data


def vrel_from_coords(x1, y1, z1, x2, y2, z2):
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2
    vrel = (dx*dx + dy*dy + dz*dz)**0.5
    return vrel


def option_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result.add_argument("setname", help="set name", type=str)

    result.add_argument("-liM", help="mass lower limit, in MSun",
                        dest="Mmin", type=float, default=10)
    result.add_argument("-lsM", help="mass upper limit, in MSun",
                        dest="Mmax", type=float, default=50)
    result.add_argument("-lia", help="minimum inner semimajor axis, in au",
                        dest="ainn_min", type=float, default=0.1)
    result.add_argument("-lsa", help="maximum inner semimajor axis, in au",
                        dest="ainn_max", type=float, default=1)
    result.add_argument("-vinf", help="velocity at infinity scale parameter, in km/s",
                        dest="sigma_vinf", type=float, default=1)
    result.add_argument("-bmax", help="maximum impact parameter, in units of ainn",
                        dest="bmax_ainn", type=float, default=2)
    result.add_argument("-Rschw", help="set collision radius as this multiple of Schwarzschild radii",
                        dest="Rschw_mult", type=float, default=10)
    result.add_argument("-seed", help="seed",
                        dest="seed", type=int, default=None)

    result = run_options(result)
    result.add_argument("-N", help="number of realizations",
                        dest="N", type=float, default=20)
    result.add_argument("-j", help="default number of cores to use in the bash script, 0 is as many as possible",
                        dest="ncores", type=int, default=8)
    result.add_argument("-RW", help="rewrite parallel file with different options",
                        dest="rewrite", type=str, default=None)
    result.add_argument("-K", help="use keplerian initial conditions",
                        dest="keplerian", type=str, default=None)
    return result


def generate_initial_condition_set(setname, N, tfin_mult, idist_mult,
                                   logging, individual, rewrite=False, keplerian=False,
                                   save_traj=False, nocoll=False, nogw=False, seed=None, ncores=8,
                                   ic_function=initial_conditions_isol, ic_arguments={}):

    input_file = setname + ".h5" if not setname.endswith(".h5") else setname
    cmd_sh_file = setname + ".sh"
    output_file = input_file.replace(".h5", output_suff)
    coll_folder = input_file.replace(".h5", coll_folder_suff)
    N = int(N)

    if rewrite is None:
        if os.path.exists(input_file):
            print("Found existing input file. Deleting it and all associated files")
            os.remove(input_file)
        if os.path.exists(cmd_sh_file): os.remove(cmd_sh_file)
        if os.path.exists(output_file): os.remove(output_file)
        if os.path.exists(coll_folder): shutil.rmtree(coll_folder)
        np.random.seed(seed)

        if keplerian is None:
            data = []
            for j in range(N):
                value_init = ic_function(**ic_arguments)
                data.append(value_init)

            datah5 = pd.DataFrame(data)
        else:
            kepler_file = keplerian + ".h5" if not keplerian.endswith(".h5") else keplerian
            datah5_kepl = pd.read_hdf(kepler_file)
            vinf = vrel_from_coords(datah5_kepl["vx_sin"], datah5_kepl["vy_sin"], datah5_kepl["vz_sin"],
                                    datah5_kepl["vx_bin"], datah5_kepl["vy_bin"], datah5_kepl["vz_bin"])

            datah5 = datah5_kepl.assign(vinf=vinf)
            datah5 = datah5.assign(b=np.fabs(datah5_kepl["bmod"]),
                                   i_out=np.fabs(datah5_kepl["Vphi_sin"]))
            datah5 = datah5.assign(ome_out=np.random.uniform(0, 2*np.pi, len(datah5_kepl)),
                                   Ome_out=np.random.uniform(0, 2*np.pi, len(datah5_kepl)))
        datah5.to_hdf(input_file, key="initial_conditions")
        N = len(datah5)

    else:
        rewrite = rewrite.replace(".h5", "")
        input_file2 = rewrite + ".h5"
        cmd_sh_file = rewrite + ".sh"
        shutil.copy(input_file, input_file2)
        input_file = input_file2
        coll_folder = input_file.replace(".h5", coll_folder_suff)
        output_file = input_file.replace(".h5", output_suff)
        datah5 = pd.read_hdf(input_file2)
        N = len(datah5)

    os.mkdir(coll_folder)

    with open(cmd_sh_file, 'w') as f:
        string = "# python3 binsingle_run.py {:s} -id {:d} -idist {:g} -tfin {:g}" \
                 "\n# Run this file with *bash* (not sh or zsh or else)\n" \
                 "\n".format(input_file, 0, idist_mult, tfin_mult)

        f.write("#!/bin/bash\n\n")
        f.write(string)
        f.write("setname={:s}\n".format(input_file.replace(".h5", "")))
        f.write("cmd='python3 binsingle_run.py'\n")
        f.write("N={:d}\n".format(N))
        f.write("tfin={:g}\n".format(tfin_mult))
        f.write("idist={:g}\n".format(idist_mult))
        f.write("outfile=$setname'_output.txt'\n")
        f.write("exmap=$setname'_exmap.txt'\n")
        f.write("ncores=${{1:-{:d}}}\n\n".format(ncores))

        f.write("logging={:d}\n".format(logging))
        f.write("individual={:d}\n".format(individual))
        f.write("save_traj={:d}\n".format(save_traj))
        f.write("nocoll={:d}\n".format(nocoll))
        f.write("nogw={:d}\n".format(nogw))

        f.write("[ $logging == 1 ] && L=' -L' || L=''\n")
        f.write("[ $individual == 1 ] && T=' -T' || T=''\n")
        f.write("[ $save_traj == 1 ] && S=' -S' || S=''\n")
        f.write("[ $nocoll == 1 ] && NC=' -nocoll' || NC=''\n")
        f.write("[ $nogw == 1 ] && NGW=' -nogw' || NGW=''\n\n")

        f.write("cmdfile=$setname'_cmdrun.txt'\n")
        f.write("[ -e $cmdfile ] && rm $cmdfile\n\n")
        f.write("for ((i = 0; i < $N; i++)); do\n\techo \"$cmd $setname.h5 -id $i -idist $idist "
                "-tfin $tfin $L$T$S$NC$NGW\" >> $cmdfile\ndone\n\n")

        f.write("[ -e $outfile ] && rm $outfile\n")
        f.write("[ -e $exmap ] && rm $exmap\n")
        f.write("parallel -a $cmdfile -j $ncores\n")
        f.write("python3 binsingle_txt_to_h5.py $outfile -isol\n".format(output_file))

    print("Written initial condition file {:s}\nRun set with {:s}".format(input_file, cmd_sh_file))


if __name__ == "__main__":
    args = option_parser().parse_args()

    ic_func = initial_conditions_isol
    ic_args = {"Mmin": args.Mmin, "Mmax": args.Mmax, "ainn_max": args.ainn_max, "ainn_min": args.ainn_min,
               "bmax_ainn": args.bmax_ainn, "sigma_vinf": args.sigma_vinf, "Rschw_mult": args.Rschw_mult}
    generate_initial_condition_set(setname=args.setname, N=args.N, tfin_mult=args.tfin_mult, idist_mult=args.idist_mult,
                                   logging=args.logging, individual=args.individual, rewrite=args.rewrite, keplerian=args.keplerian,
                                   save_traj=args.save_traj, nocoll=args.nocoll, nogw=args.nogw, seed=args.seed, ncores=args.ncores,
                                   ic_function=ic_func, ic_arguments=ic_args)

