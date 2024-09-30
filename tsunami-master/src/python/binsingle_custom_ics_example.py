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

Example usage:
    * python binsingle_initial_conditions.py testset -N 100
    * bash testset.sh

Creates a set named testset, with 100 binary-single encounters, using the default initial condition parameters
Check the varius input parameters in the argument options with -h
"""


import tsunami
import numpy as np
from binsingle_initial_conditions import generate_initial_condition_set
import argparse

KU = tsunami.KeplerUtils()

output_suff = "_output.txt"
coll_folder_suff = "_coll"


def initial_conditions_chaos(b_ainn=1, ainn_min=0.1, ainn_max=1, vinf=1.0):
    m1 = 17.5
    m2 = 15
    m3 = 12.5

    R1 = 0.0
    R2 = 0.0
    R3 = 0.0

    a_inn = 5
    b = a_inn * b_ainn
    e_inn = 0.0
    ome_inn = 0.0

    v_to_kms = 2 * np.pi * KU.au2km / KU.yr2sec
    vinf_nb = vinf / v_to_kms

    M_inn = np.random.uniform(0, 2 * np.pi)
    nu_inn = KU.M_to_nu(M_inn, e_inn)

    i_out = np.random.uniform(-1, 1)
    i_out = np.arccos(i_out)
    Ome_out = 0.0
    ome_out = 0.0

    return m1, m2, m3, R1, R2, R3, a_inn, e_inn, ome_inn, nu_inn, vinf_nb, b, i_out, ome_out, Ome_out


def option_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    result.add_argument("setname", help="set name",
                        type=str)

    result.add_argument("-N", help="number of realizations",
                        dest="N", type=float, default=20)
    result.add_argument("-j", help="default number of cores to use in the bash script, 0 is as many as possible",
                        dest="ncores", type=int, default=8)
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
                        dest="bmax_ainn", type=float, default=1)
    result.add_argument("-seed", help="seed",
                        dest="seed", type=int, default=None)

    result.add_argument("-tfin", help="final time, in units of time_to_encounter",
                        dest="tfin_mult", type=float, default=1e3)
    result.add_argument("-idist", help="initial distance, in units of binary semimajor axis",
                        dest="idist_mult", type=float, default=150)
    result.add_argument("-Rschw", help="collision radius, in units of Schwarzschild radii",
                        dest="Rschw_mult", type=float, default=10)
    result.add_argument("-S", help="save trajectories",
                        dest="save_traj", action="store_true")
    result.add_argument("-L", help="save logfile",
                        dest="logging", action="store_true")
    result.add_argument("-I", help="individual save",
                        dest="individual", action="store_true")
    result.add_argument("-nogw", help="do not use gw",
                        dest="nogw", action="store_true")
    result.add_argument("-nocoll", help="do not use collisions",
                        dest="nocoll", action="store_true")

    result.add_argument("-RW", help="rewrite parallel file with different options",
                        dest="rewrite", type=str, default=None)
    result.add_argument("-K", help="use keplerian initial conditions",
                        dest="keplerian", type=str, default=None)
    return result


if __name__ == "__main__":
    args = option_parser().parse_args()

    ic_func = initial_conditions_chaos
    ic_args = {"b_ainn": 1.0, "ainn_min": 0.1, "ainn_max": 1, "vinf": 1.0}
    generate_initial_condition_set(setname=args.setname, N=args.N, tfin_mult=args.tfin_mult, idist_mult=args.idist_mult,
                                   logging=args.logging, individual=args.individual, rewrite=args.rewrite, keplerian=args.keplerian,
                                   save_traj=args.save_traj, nocoll=args.nocoll, nogw=args.nogw, seed=args.seed, ncores=args.ncores,
                                   ic_function=ic_func, ic_arguments=ic_args)
