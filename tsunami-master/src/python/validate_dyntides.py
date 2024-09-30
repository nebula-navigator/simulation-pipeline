import tsunami
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

keplutils = tsunami.KeplerUtils()


def setup_binary(a=-1.0, e=0.99, m1=50.0, m2=50.0, nu=0.0, R1=40, R2=40, d0=None, rp=None):
    i = ome = Ome = 0.0 # rad
    pos_vel2 = np.array([0.,0.,0., 0.,0.,0.])

    if d0 is not None and rp is not None:
        rp = rp * keplutils.RSun2au
        e = 1 - rp / a
        p = a * (1 - e * e)
        nu = p / d0 - 1
        nu = nu / e
        nu = - np.arccos(nu)
    else: rp = a*(1-e)

    print("\nGenerating binary")
    print("m1 = {:g} MSun".format(m1))
    print("m2 = {:g} MSun".format(m2))
    print("a = {:g} au".format(a))
    print("e = {:g}".format(e))
    print("rp = {:g} au".format(rp))
    if a < 0:
        vinf = (-(m1+m2) / a)**0.5
        vinf = vinf * 2*np.pi * keplutils.au2km / keplutils.yr2sec
        print("vinf = {:g} km/s".format(vinf))
    print("nu = {:g} deg".format(np.degrees(nu)))

    R = np.array([R1, R2], dtype=np.float64)

    R = R * keplutils.RSun2au
    print("Radii (au)", R)
    pos_vel1 = keplutils.kepl_to_cart(pos_vel2, m1, m2, a, e, i, ome, Ome, nu)

    m = np.array([m1, m2])
    p = np.array([pos_vel1[:3], pos_vel2[:3]])
    v = np.array([pos_vel1[3:], pos_vel2[3:]])

    return m, p, v, R


def run_binary(m, p, v, R, poly1=1.5, poly2=1.5, ft=500, dt=1, datname="validate_dyntides"):
    code = tsunami.Tsunami()  # units: G=1, L=au, M=MSun
    # Use dynamical tides
    code.wDynTides = True


    st = np.array([0, 0])
    kaps = np.array([0.0, 1.0])
    taus = np.array([1e-10, 1e-10])
    poly = np.array([poly1, poly2])

    print(p, v, m, R)

    code.add_particle_set(p, v, m, R, st)
    code.initialize_tides(kaps, taus, poly)

    print("tfin = {:g} yr".format(ft))
    print("dt = {:g} yr".format(dt))

    tfin = ft / code.Tscale
    t = 0.0
    dt = dt / code.Tscale

    t_l, e_l, a_l = [], [], []
    print("\nEvolving binary")

    while t < tfin:
        code.evolve_system(t)
        t = code.ctime  # synchronize to internal time

        code.sync_position_and_velocities(p, v)

        if code.CollLogger.collflag:
            print("\nCollision at time {:g} yr".format(t * code.Tscale))
            break

        dp = p[0]-p[1]
        dv = v[0]-v[1]
        orb = keplutils.cart_to_kepl(dp, dv, m[0], m[1])

        print("Distance / Radius :", (dp*dp).sum()**0.5 / R.sum())

        t_l.append(t * code.Tscale)
        a_l.append(orb[0])
        e_l.append(orb[1])

        t += dt

        print("Time: {:4.4f} yr, t/tfin {:2.1%}. Delta Err = {:1.4e}".format(code.ctime * code.Tscale, code.ctime/tfin, code.deltaE), end="\r")

    print("Saving evolution to", datname + ".txt\033[K")
    np.savetxt(datname + ".txt", np.vstack((t_l, a_l, e_l)).T)


def plot_evolution(datname="validate_dyntides"):
    data = np.loadtxt(datname + ".txt")
    t_l, a_l, e_l = data.T

    f, ax = plt.subplots(2, sharex=True, figsize=(7, 5))

    ax[0].plot(t_l, a_l, ls="-", lw=3)
    ax[0].set_ylabel("semimajor axis [au]")

    ax[1].plot(t_l, np.array(e_l), label="", ls="-", lw=3)
    ax[1].set_ylabel("eccentricity")
    ax[1].set_xlabel("time [yr]")
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
    result.add_argument("-a", dest="a", default=-1.0,
                        type=float, help="semimajor axis in au")
    exclude = result.add_mutually_exclusive_group()
    exclude.add_argument("-e", dest="e", default=None,
                        type=float, help="eccentricity")
    exclude.add_argument("-rp", dest="rp", default=500,
                         type=float, help="pericenter distance in RSun")
    result.add_argument("-m1", dest="m1", default=50.0,
                        type=float, help="mass object 1 in MSun")
    result.add_argument("-m2", dest="m2", default=50.0,
                        type=float, help="mass object 2 in MSun")
    result.add_argument("-R1", dest="R1", default=40.0,
                        type=float, help="radius object 1 in RSun")
    result.add_argument("-R2", dest="R2", default=40.0,
                        type=float, help="radius object 2 in RSun")
    result.add_argument("-poly1", dest="poly1", default=1.5,
                        type=float, help="polytropic index object 1")
    result.add_argument("-poly2", dest="poly2", default=1.5,
                        type=float, help="polytropic index object 2")
    exclude = result.add_mutually_exclusive_group()
    exclude.add_argument("-nu", dest="nu", default=None,
                        type=float, help="initial binary true anomaly, in radians")
    exclude.add_argument("-d", dest="d0", default=1,
                        type=float, help="initial distance of hyperbolic orbit, in au")
    result.add_argument("-ft", dest="ft", default=160,
                        type=float, help="final time in yr")
    result.add_argument("-dt", dest="dt", default=1,
                        type=float, help="output timestep in yr")
    result.add_argument("-dat", dest="datname", default="validate_dyntides",
                        type=str, help="text file with evolution data")
    return result


if __name__ == "__main__":
    args = new_argument_parser().parse_args()
    print("Validating post-Newtonians, may take a few minutes")

    if not (os.path.isfile(args.datname+".txt")):
        m, p, v, R = setup_binary(args.a, args.e, args.m1, args.m2, args.nu, args.R1, args.R2, args.d0, args.rp)
        run_binary(m, p, v, R, args.poly1, args.poly2, args.ft, args.dt, args.datname)
    else:
        print("Found {:s}, skipping simulation and plotting data from file".format(args.datname))
    plot_evolution(args.datname)
