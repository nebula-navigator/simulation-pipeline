import numpy as np
import argparse
import signal
import h5py
import tsunami
import os


def period(m, a):
    P = 2 * np.pi * (a * a * a / m) ** 0.5
    return P


def semi(mx, my, P):
    mu = mx + my
    a = (((P / (2.0 * np.pi)) ** 2) * mu) ** (1.0 / 3.0)
    return a


def common_option_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    result.add_argument("-m1", help="m1, MSun",
                        dest="m1", type=float, default=0.01)
    result.add_argument("-m2", help="m2, MSun",
                        dest="m2", type=float, default=1)
    result.add_argument("-m3", help="m3, MSun",
                        dest="m3", type=float, default=0.04)
    result.add_argument("-R1", help="R1, RSun",
                        dest="R1", type=float, default=1e-3)
    result.add_argument("-R2", help="R2, RSun",
                        dest="R2", type=float, default=1.0)
    result.add_argument("-R3", help="R3, RSun",
                        dest="R3", type=float, default=0.04)

    exclude1 = result.add_mutually_exclusive_group()
    exclude1.add_argument("-P1", help="P1, yr",
                          dest="P1", type=float, default=None)
    exclude1.add_argument("-a1", help="a1, au",
                          dest="a1", type=float, default=6)
    result.add_argument("-e1", help="e1",
                        dest="e1", type=float, default=0.001)

    exclude2 = result.add_mutually_exclusive_group()
    exclude2.add_argument("-P2", help="P2, yr",
                          dest="P2", type=float, default=None)
    exclude2.add_argument("-a2", help="a2, au",
                          dest="a2", type=float, default=100)
    result.add_argument("-e2", help="e2",
                        dest="e2", type=float, default=0.6)

    result.add_argument("-imut", help="mutual inclination, degrees",
                        dest="i_mut", type=float, default=65)
    result.add_argument("-ome1", help="ome1, degrees",
                        dest="ome1", type=float, default=45)
    result.add_argument("-ome2", help="ome2, degrees",
                        dest="ome2", type=float, default=0)
    result.add_argument("-Ome2", help="Ome2, degrees",
                        dest="Ome2", type=float, default=0.0)

    result.add_argument("-ft", help="final time, yr",
                        dest="ftime", type=float, default=2.5e7)
    result.add_argument("-dt_out", help="output timestep, yr",
                        dest="dt_out", type=float, default=1e4)

    result.add_argument("-k1", help="apsidal motion constant of particle 1",
                        dest="k1", type=float, default=0.014)
    result.add_argument("-k2", help="apsidal motion constant of particle 2",
                        dest="k2", type=float, default=0.3)
    result.add_argument("-k3", help="apsidal motion constant of particle 3",
                        dest="k3", type=float, default=0.3)
    result.add_argument("-tau1", help="time-lag of particle 1, seconds",
                        dest="tau1", type=float, default=100e4)
    result.add_argument("-tau2", help="time-lag of particle 2, seconds",
                        dest="tau2", type=float, default=0.66e3)
    result.add_argument("-tau3", help="time-lag of particle 3, seconds",
                        dest="tau3", type=float, default=0.66e3)
    result.add_argument("-gyr1", help="gyration radius of particle 1",
                        dest="gyr1", type=float, default=0.28)
    result.add_argument("-gyr2", help="gyration radius of particle 2",
                        dest="gyr2", type=float, default=0.28)
    result.add_argument("-gyr3", help="gyration radius of particle 3",
                        dest="gyr3", type=float, default=0.28)

    result.add_argument("-Prot1", help="rotational period of particle 1, yr",
                        dest="Prot1", type=float, default=0.3)
    result.add_argument("-Prot2", help="rotational period of particle 2, yr",
                        dest="Prot2", type=float, default=0.3)
    result.add_argument("-Prot3", help="rotational period of particle 3, yr",
                        dest="Prot3", type=float, default=0.3)

    result.add_argument("-obl1", help="obliquity of particle 1, degrees",
                        dest="obl1", type=float, default=30)
    result.add_argument("-obl2", help="obliquity of particle 2, degrees",
                        dest="obl2", type=float, default=120)
    result.add_argument("-obl3", help="obliquity of particle 3, degrees",
                        dest="obl3", type=float, default=0.0)

    result.add_argument("-nocol", help="no collision",
                        dest="nocol", action="store_true")
    result.add_argument("-tide", help="use tides",
                        dest="tide", action="store_true")
    result.add_argument("-gr", help="use post-Newtonian corrections",
                        dest="gr", action="store_true")
    result.add_argument("-nopn1", help="turn off PN1",
                        dest="nopn1", action="store_true")
    result.add_argument("-nopn2", help="turn off PN2",
                        dest="nopn2", action="store_true")
    result.add_argument("-nopn25", help="turn off PN2.5",
                        dest="nopn25", action="store_true")
    result.add_argument("-nopn3", help="turn off PN3",
                        dest="nopn3", action="store_true")
    result.add_argument("-nopn35", help="turn off PN3.5",
                        dest="nopn35", action="store_true")

    result.add_argument("--mute", help="mute output (do not save outputs)",
                        dest="mute_output", action="store_true")
    result.add_argument("--R", help="restart simulation, erasing previous evolution",
                        dest="restart", action="store_true")

    result.add_argument("-task", "--t", dest="task", default=['show'], nargs='+', type=str,
                        help="task for output values")
    result.add_argument("filename", type=str, default="test_triple",
                        help="output hdf5 file", nargs='?')
    return result


def okinami_option_parser():
    result = common_option_parser()

    result.add_argument("-stabcrit", help="stability criterion",
                        dest="stabcrit", type=int, default=0)
    result.add_argument("-qsec", help="stop if quasi-secular evolution is detected",
                        dest="qsec", action="store_true")
    return result


def tsunami_option_parser():
    result = common_option_parser()

    result.add_argument("-nu1", help="nu1, degrees",
                        dest="nu1", type=float, default=0.0)
    result.add_argument("-nu2", help="nu2, degrees",
                        dest="nu2", type=float, default=135)
    result.add_argument("--invariant", "--inv", help="set up system in the invariant plane",
                        action="store_true", dest="invariant")
    return result


def process_args(args):
    """
    Here we convert everything to N-body units / radians
    :param args:
    :return:
    """

    KU = tsunami.KeplerUtils()

    ome1, ome2, Ome2 = map(lambda x: np.random.uniform(0, 2 * np.pi) if x is None else np.radians(x),
                           [args.ome1, args.ome2, args.Ome2])

    if args.P1 is not None:
        a1 = semi(args.m1, args.m2, args.P1 * KU.Tscale)
    else:
        a1 = args.a1

    if args.P2 is not None:
        a2 = semi(args.m1, args.m3, args.P2 * KU.Tscale)
    else:
        a2 = args.a2

    i_mut = np.radians(args.i_mut)

    m1, m2, m3 = args.m1, args.m2, args.m3
    e1, e2 = args.e1, args.e2
    R1, R2, R3 = args.R1 * KU.RSun2au, args.R2 * KU.RSun2au, args.R3 * KU.RSun2au
    k1, k2, k3 = args.k1, args.k2, args.k3
    tau1, tau2, tau3 = args.tau1 / KU.yr2sec / KU.Tscale, args.tau2 / KU.yr2sec / KU.Tscale, args.tau3 / KU.yr2sec / KU.Tscale

    spin1, spin2, spin3 = map(lambda x: 2 * np.pi / x * KU.Tscale, [args.Prot1, args.Prot2, args.Prot3])
    obl1, obl2, obl3 = map(lambda x:  np.radians(x), [args.obl1, args.obl2, args.obl3])

    return m1, m2, m3, R1, R2, R3, a1, e1, a2, e2, i_mut, ome1, ome2, Ome2, k1, k2, k3, tau1, tau2, tau3, \
        spin1, spin2, spin3, obl1, obl2, obl3


def update_hdf5(hdfile, tag, data, t=np.double):
    if tag in hdfile:
        hdfile[tag][...] = data
    else:
        hdfile[tag] = np.asarray(data, dtype=t)


class TripleSystem():
    hdf_ic_nbody = "ic_nbody"
    hdf_ic_secular = "ic_secular"
    hdf_ic_kepler = "ic_kepler"
    hdf_time_data = "time_data"

    def __init__(self, m1, m2, m3, a1, a2, e1, e2, i_mut, ome1, ome2, Ome2, R1, R2, R3,
                 k1=0.0, k2=0.0, k3=0.0, tau1=0.0, tau2=0.0, tau3=0.0, gyr1=0.0, gyr2=0.0, gyr3=0.0,
                 spin1=0.0, spin2=0.0, spin3=0.0, obl1=0.0, obl2=0.0, obl3=0.0, azim1=0.0, azim2=0.0, azim3=0.0,
                 nopn1=False, nopn2=False, nopn25=False, nopn3=False, nopn35=False, nocol=True, tide=False,
                 gr=False, filename="test_triple", mute_output=False, restart=False):
        """
        All units are in N-body units (G=1, L=au, M=MSun, T=yr/(2pi) already
        ***Obliquity here is with respect the x-y plane, not the plane of the binary***
        ***For invariant_plane=False (default) the inner binary is in the x-y plane anyway***
        :param m1:
        :param m2:
        :param m3:
        :param a1:
        :param a2:
        :param e1:
        :param e2:
        :param i_mut:
        :param ome1:
        :param ome2:
        :param Ome2:
        :param R1:
        :param R2:
        :param R3:
        :param k1:
        :param k2:
        :param k3:
        :param tau1:
        :param tau2:
        :param tau3:
        :param gyr1:
        :param gyr2:
        :param gyr3:
        :param spin1:
        :param spin2:
        :param spin3:
        :param obl1:
        :param obl2:
        :param obl3:
        :param azim1:
        :param azim2:
        :param azim3:
        :param nocol:
        :param tide:
        :param gr:
        :param filename:
        :param mute_output:
        :param restart:
        """

        self.a1_0 = self.a1 = a1
        self.a2_0 = self.a2 = a2

        self.m1, self.m2, self.m3 = m1, m2, m3
        self.e1_0 = self.e1 = e1
        self.e2_0 = self.e2 = e2

        self.ome1_0 = self.ome1 = ome1
        self.ome2_0 = self.ome2 = ome2

        self.Ome2_0 = self.Ome2 = Ome2
        self.i_mut_0 = self.i_mut = i_mut

        self.R1, self.R2, self.R3 = R1, R2, R3
        self.k1, self.k2, self.k3 = k1, k2, k3

        self.KU = tsunami.KeplerUtils()
        self.time_unit = 1.0 / self.KU.Tscale
        self.tau1, self.tau2, self.tau3 = tau1, tau2, tau3
        self.gyr1, self.gyr2, self.gyr3 = gyr1, gyr2, gyr3

        self.spin1, self.spin2, self.spin3 = spin1, spin2, spin3
        self.obl1, self.obl2, self.obl3 = obl1, obl2, obl3
        self.azim1, self.azim2, self.azim3 = azim1, azim2, azim3

        self.nocol = nocol
        self.tide = tide
        self.gr = gr
        self.mute_output = mute_output
        self.restart = restart
        self.nopn1, self.nopn2, self.nopn25, self.nopn3, self.nopn35 = nopn1, nopn2, nopn25, nopn3, nopn35

        signal.signal(signal.SIGINT, self.exit_sig)
        signal.signal(signal.SIGTERM, self.exit_sig)
        self.killed = False
        self.skip_run = False
        self.secularcode = self.directcode = None
        self.skip_file = False

        self.outname = filename
        if not self.outname.endswith(".hdf5"): self.outname = self.outname + ".hdf5"

    def info_system(self):
        if not self.mute_output:
            print("Filename: {:s}".format(os.path.basename(self.outname)))
        else:
            print("No filename in output")
        print("Path: {:s}".format(os.path.dirname(self.outname)))

        print("m1, m2, m3 = {:g} MSun, {:g} MSun, {:g} MSun".format(self.m1, self.m2, self.m3))
        print("a1, e1 = {:g} au, {:g}".format(self.a1_0, self.e1_0))
        print("a2, e2 = {:g} au, {:g}".format(self.a2_0, self.e2_0))
        print("i_mut = {:g} deg".format(np.degrees(self.i_mut)))

        self.P1 = period(self.m1 + self.m2, self.a1_0)
        self.P2 = period(self.m1 + self.m2 + self.m3, self.a2_0)
        print("P1 = {:g} yr".format(self.P1 / self.time_unit))
        print("P2 = {:g} yr".format(self.P2 / self.time_unit))

        print("R1, R2, R3 (au):", self.R1, self.R2, self.R3)
        print("k1, k2, k3:", self.k1, self.k2, self.k3)
        print("tau1, tau2, tau3 (nbody):", self.tau1, self.tau2, self.tau3)

        if self.directcode is not None and self.directcode.Conf.wSpins:
            if self.spin1 != 0: print("Prot1 = {:g} yr".format(2 * np.pi / self.spin1 * self.KU.Tscale))
            print("spin1 = {:g} 2pi/yr".format(self.spin1 / self.KU.Tscale))
            if self.spin2 != 0: print("Prot2 = {:g} yr".format(2 * np.pi / self.spin2 * self.KU.Tscale))
            print("spin2 = {:g} 2pi/yr".format(self.spin2 / self.KU.Tscale))
            if self.spin3 != 0: print("Prot3 = {:g} yr".format(2 * np.pi / self.spin3 * self.KU.Tscale))
            print("spin3 = {:g} 2pi/yr".format(self.spin3 / self.KU.Tscale))

        print("GR corrections:", self.gr)
        print("Tides:", self.tide)
        print("Collisions:", not self.nocol)

    def setup_spin_vectors(self, spin1, spin2, spin3, obl1, obl2, obl3,
                           azim1=0.0, azim2=0.0, azim3=0.0):

        spin1_vec = spin1 * np.array([np.sin(obl1) * np.sin(azim1),
                                           -np.sin(obl1) * np.cos(azim1),
                                           np.cos(obl1)])

        spin2_vec = spin2 * np.array([np.sin(obl2) * np.sin(azim2),
                                           -np.sin(obl2) * np.cos(azim2),
                                           np.cos(obl2)])

        spin3_vec = spin3 * np.array([np.sin(obl3) * np.sin(azim3),
                                           -np.sin(obl3) * np.cos(azim3),
                                           np.cos(obl3)])
        s = np.array([spin1_vec, spin2_vec, spin3_vec], dtype=np.float64)
        return s

    def exit_sig(self, signum, frame):
            print("Signal {:d} received, exiting".format(signum))
            self.killed = True

    def save_record_par(self, hdfile,
                        recordpar, npar):
        if self.hdf_time_data + "/recordpar" in hdfile:
            npar_old = hdfile[self.hdf_time_data].attrs["npar"]
            if npar_old != npar:
                raise ValueError("Old npar ({:d}) different than new one {:d}".format(npar_old, npar))

            del hdfile[self.hdf_time_data].attrs["npar"]
            del hdfile[self.hdf_time_data + "/recordpar"]

        hdfile[self.hdf_time_data + "/recordpar"] = recordpar.astype('float32')

        hdfile[self.hdf_time_data].attrs["npar"] = npar

    def load_record_par(self, hdfile):
        recordpar = hdfile[self.hdf_time_data + "/recordpar"][()] if self.hdf_time_data + "/recordpar" in hdfile else []
        return recordpar

    def save_kepler_ic(self, hdfile,
                       m1, m2, m3, a1, a2, e1, e2, i_mut,
                       ome1, ome2, Ome2, R1, R2, R3, k1, k2, k3, tau1, tau2, tau3, nu1, nu2, gyr1, gyr2, gyr3,
                       spin1, spin2, spin3, obl1, obl2, obl3, azim1, azim2, azim3):
        update_hdf5(hdfile, self.hdf_ic_kepler + "/m1", m1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/m2", m2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/m3", m3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/a1", a1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/a2", a2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/e1", e1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/e2", e2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/R1", R1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/R2", R2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/R3", R3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/ome1", ome1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/ome2", ome2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/Ome2", Ome2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/i_mut", i_mut)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/k1", k1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/k2", k2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/k3", k3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/tau1", tau1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/tau2", tau2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/tau3", tau3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/gyr1", gyr1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/gyr2", gyr2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/gyr3", gyr3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/spin1", spin1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/spin2", spin2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/spin3", spin3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/obl1", obl1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/obl2", obl2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/obl3", obl3)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/azim1", azim1)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/azim2", azim2)
        update_hdf5(hdfile, self.hdf_ic_kepler + "/azim3", azim3)
        if nu1 is not None: update_hdf5(hdfile, self.hdf_ic_kepler + "/nu1", nu1)
        if nu2 is not None: update_hdf5(hdfile, self.hdf_ic_kepler + "/nu2", nu2)

    def load_kepler_ic(self, hdfile):
        m1 = hdfile[self.hdf_ic_kepler + "/m1"][()]
        m2 = hdfile[self.hdf_ic_kepler + "/m2"][()]
        m3 = hdfile[self.hdf_ic_kepler + "/m3"][()]
        a1 = hdfile[self.hdf_ic_kepler + "/a1"][()]
        a2 = hdfile[self.hdf_ic_kepler + "/a2"][()]
        e1 = hdfile[self.hdf_ic_kepler + "/e1"][()]
        e2 = hdfile[self.hdf_ic_kepler + "/e2"][()]
        R1 = hdfile[self.hdf_ic_kepler + "/R1"][()]
        R2 = hdfile[self.hdf_ic_kepler + "/R2"][()]
        R3 = hdfile[self.hdf_ic_kepler + "/R3"][()]
        ome1 = hdfile[self.hdf_ic_kepler + "/ome1"][()]
        ome2 = hdfile[self.hdf_ic_kepler + "/ome2"][()]
        Ome2 = hdfile[self.hdf_ic_kepler + "/Ome2"][()]
        i_mut = hdfile[self.hdf_ic_kepler + "/i_mut"][()]
        k1 = hdfile[self.hdf_ic_kepler + "/k1"][()]
        k2 = hdfile[self.hdf_ic_kepler + "/k2"][()]
        k3 = hdfile[self.hdf_ic_kepler + "/k3"][()]
        tau1 = hdfile[self.hdf_ic_kepler + "/tau1"][()]
        tau2 = hdfile[self.hdf_ic_kepler + "/tau2"][()]
        tau3 = hdfile[self.hdf_ic_kepler + "/tau3"][()]
        gyr1 = hdfile[self.hdf_ic_kepler + "/gyr1"][()]
        gyr2 = hdfile[self.hdf_ic_kepler + "/gyr2"][()]
        gyr3 = hdfile[self.hdf_ic_kepler + "/gyr3"][()]
        spin1 = hdfile[self.hdf_ic_kepler + "/spin1"][()]
        spin2 = hdfile[self.hdf_ic_kepler + "/spin2"][()]
        spin3 = hdfile[self.hdf_ic_kepler + "/spin3"][()]
        obl1 = hdfile[self.hdf_ic_kepler + "/obl1"][()]
        obl2 = hdfile[self.hdf_ic_kepler + "/obl2"][()]
        obl3 = hdfile[self.hdf_ic_kepler + "/obl3"][()]
        azim1 = hdfile[self.hdf_ic_kepler + "/azim1"][()]
        azim2 = hdfile[self.hdf_ic_kepler + "/azim2"][()]
        azim3 = hdfile[self.hdf_ic_kepler + "/azim3"][()]
        nu1 = None if not self.hdf_ic_kepler + "/nu1" in hdfile else hdfile[self.hdf_ic_kepler + "/nu1"][()]
        nu2 = None if not self.hdf_ic_kepler + "/nu2" in hdfile else hdfile[self.hdf_ic_kepler + "/nu2"][()]
        return m1, m2, m3, a1, a2, e1, e2, i_mut, ome1, ome2, Ome2, R1, R2, R3, k1, k2, k3, tau1, tau2, tau3, nu1, nu2, \
            gyr1, gyr2, gyr3, spin1, spin2, spin3, obl1, obl2, obl3, azim1, azim2, azim3

    def plot(self, task, wSpins=False):
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ModuleNotFoundError:
            print("Sorry, you don't have seaborn or matplotlib.\nAborting plotting")
            return

        sns.set(font_scale=1.33)
        sns.set_style("ticks")

        hfont = {'fontname': 'DejaVu Serif', 'fontsize': 16, 'labelpad': 10}

        recordpar = np.asarray(self.recordpar)
        recordpar = recordpar.reshape((-1, self.npar)).T

        time = recordpar[0]
        e1 = recordpar[1]
        e2 = recordpar[2]
        a1 = recordpar[3]
        a2 = recordpar[4]
        i_mut = np.degrees(recordpar[5])

        period_from_a1 = period(self.m1 + self.m2, a1) / self.time_unit
        period_from_a2 = period(self.m1 + self.m2 + self.m3, a2) / self.time_unit

        nplots = 3
        figsize = np.array((12, 9))
        if wSpins:
            nplots+=1
            hfont["fontsize"] = 14
            figsize[1]+=3

        f, ax = plt.subplots(nplots, figsize=(12, 9), sharex=True, gridspec_kw=dict(hspace=0),  tight_layout=True)
        ax[0].plot(time, e1, label="$e_1$", lw=2.5, alpha=0.66)
        ax[0].plot(time, e2, label="$e_2$", lw=2.5, alpha=0.66)
        ax[0].set_ylabel("eccentricity", **hfont)
        ax[0].tick_params(labelsize=hfont["fontsize"])
        ax[0].legend()
        ax[1].plot(time, i_mut, lw=2.5, label="$i_{\\rm mut}$", c="black")
        ax[1].set_ylabel("inclination [deg]", **hfont)
        ax[1].set_xlabel("time [yr]", **hfont)
        ax[1].tick_params(labelsize=hfont["fontsize"])
        ax[1].legend()

        ax[2].plot(time, period_from_a1, lw=2.5, alpha=0.66)
        ax[2].plot(time, period_from_a2, lw=2.5, alpha=0.66)
        ax[2].set_yscale("log")

        ax[2].set_ylabel("period [yr]", **hfont)
        ax[2].set_xlabel("time [yr]", **hfont)
        ax[2].tick_params(labelsize=hfont["fontsize"])

        if wSpins:
            spin1 = recordpar[8]
            spin2 = recordpar[9]
            obl1 = np.degrees(recordpar[10])
            obl2 = np.degrees(recordpar[11])
            n1 = 2*np.pi / period_from_a1

            ax[1].plot(time, obl1, lw=2.5, label="$\\theta_1$", alpha=0.66)
            ax[1].plot(time, obl2, lw=2.5, label="$\\theta_2$", alpha=0.66)
            ax[1].legend()
            ax[3].plot(time, spin1, lw=2.5, label="$S_1$", alpha=0.66)
            ax[3].plot(time, spin2, lw=2.5, label="$S_2$", alpha=0.66)
            ax[3].plot(time, n1, lw=2.5, label="$n_1$", c="black")
            ax[3].set_ylabel("spin [rad/yr]", **hfont)
            ax[3].tick_params(labelsize=hfont["fontsize"])
            ax[3].set_yscale("log")
            ax[3].set_xlabel("time [yr]", **hfont)
            ax[3].legend()

        for axx in ax:
            axx.margins(x=1e-2)
            axx.grid()

        if "show" in task: plt.show()
        if "png" in task: plt.savefig(self.outname + ".png")
        return f
