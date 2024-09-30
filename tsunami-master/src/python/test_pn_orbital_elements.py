import tsunami
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

keplutils = tsunami.KeplerUtils()


def setup_binary(a=0.1, e=0.99, m1=50.0, m2=50.0, nu=0.0, pn1=True, pn2=True):
    i = ome = Ome = 0.0 # rad
    pos_vel2 = np.array([0.,0.,0., 0.,0.,0.])

    print("\nGenerating binary")

    print("kepl_to_cart corrections:\n  PN1 = {}\n  PN2 = {}".format(pn1, pn2))

    pos_vel1 = keplutils.kepl_to_cart(pos_vel2, m1, m2, a, e, i, ome, Ome, nu, pn1=pn1, pn2=pn2)

    m = np.array([m1, m2])
    p = np.array([pos_vel1[:3], pos_vel2[:3]])
    v = np.array([pos_vel1[3:], pos_vel2[3:]])

    return m, p, v


def run_binary(m, p, v, ft=500, dt=1, datname="validate_pn", pn1=True, pn2=True):
    code = tsunami.Tsunami()  # units: G=1, L=au, M=MSun
    # Use 1PN + 2PN + 2.5PN
    code.Conf.wPNs = True

    R = 10* 2*m / code.speed_of_light**2
    st = np.array([0, 0])

    code.add_particle_set(p, v, m, R, st)

    tfin = ft / code.Tscale
    t = 0.0
    dt = dt / code.Tscale

    t_l, e_l, a_l = [], [], []
    print("\nEvolving binary")
    print("cart_to_kepl corrections:\n  PN1 = {}\n  PN2 = {}".format(pn1, pn2))

    while t < tfin:
        t += dt
        code.evolve_system(t)
        t = code.time  # synchronize to internal time

        code.sync_internal_state(p, v)

        if code.check_collision():
            print("Collision at time {:g} yr".format(t * code.Tscale))
            break

        dp = p[0]-p[1]
        dv = v[0]-v[1]
        orb = keplutils.cart_to_kepl(dp, dv, m[0], m[1], pn1=pn1, pn2=pn2)

        t_l.append(t * code.Tscale)
        a_l.append(orb[0])
        e_l.append(orb[1])

        print("Time: {:4.4f} yr, t/tfin {:2.1%}. Delta Err = {:1.4e}".format(code.time * code.Tscale, code.time/tfin, code.deltaE), end="\r")

    print("Saving evolution to", datname + ".txt\033[K")
    np.savetxt(datname + ".txt", np.vstack((t_l, a_l, e_l)).T)


def plot_evolution(datname="validate_pn"):
    data = np.loadtxt(datname + ".txt")
    t_l, a_l, e_l = data.T

    data = np.loadtxt(datname + "_pncorr.txt")
    tpn_l, apn_l, epn_l = data.T

    f, ax = plt.subplots(2, sharex=True, figsize=(7, 5))

    ax[0].plot(tpn_l, apn_l, lw=3)
    ax[0].plot(t_l, a_l, ls=":", lw=3)
    ax[0].set_ylabel("semimajor axis [au]")

    ax[1].plot(tpn_l, 1-np.array(epn_l), label="with PNs corrections", lw=3)
    ax[1].plot(t_l, 1-np.array(e_l), label="without PNs corrections", ls=":", lw=3)
    ax[1].set_ylabel("1-eccentricity")
    ax[1].set_xlabel("time [yr]")
    ax[1].set_yscale("log")
    ax[1].legend()

    plt.tight_layout()
    for axx in ax: axx.grid(True)

    f.subplots_adjust(hspace=0)

    figname = datname + ".png"
    print("Saving figure to", figname)
    plt.savefig(figname)
    plt.show()


def new_argument_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result.add_argument("-a", dest="a", default=0.1,
                        type=float, help="semimajor axis in au")
    result.add_argument("-e", dest="e", default=0.99,
                        type=float, help="eccentricity")
    result.add_argument("-m1", dest="m1", default=50.0,
                        type=float, help="mass object 1 in MSun")
    result.add_argument("-m2", dest="m2", default=50.0,
                        type=float, help="mass object 2 in MSun")
    result.add_argument("-nu", dest="nu", default=0.0,
                        type=float, help="initial binary true anomaly, in radians")
    result.add_argument("-ft", dest="ft", default=160,
                        type=float, help="final time in yr")
    result.add_argument("-dt", dest="dt", default=1,
                        type=float, help="output timestep in yr")
    result.add_argument("-dat", dest="datname", default="validate_pn",
                        type=str, help="text file with evolution data")
    return result


if __name__ == "__main__":
    args = new_argument_parser().parse_args()
    print("Validating post-Newtonians, may take a few minutes")

    if not (os.path.isfile(args.datname+".txt") and os.path.isfile(args.datname+"_pncorr.txt")):
        print("Binary:\na = {:g} au\ne = {:g}\nm1 = {:g} MSun\nm2 = {:g} MSun\nnu = {:g}".format(args.a, args.e, args.m1, args.m2, args.nu))
        m, p, v = setup_binary(args.a, args.e, args.m1, args.m2, args.nu, pn1=True, pn2=True)
        run_binary(m, p, v, args.ft, args.dt, args.datname+"_pncorr", pn1=True, pn2=True)

        m, p, v = setup_binary(args.a, args.e, args.m1, args.m2, args.nu, pn1=False, pn2=False)
        run_binary(m, p, v, args.ft, args.dt, args.datname, pn1=False, pn2=False)
    else:
        print("Found {:s}, skipping simulation and plotting data from file".format(args.datname))
    plot_evolution(args.datname)
