import numpy as np
import warnings
import h5py
import os
import tsunami
from tsunami import stateid
from triple_generator_common import okinami_option_parser, process_args, update_hdf5, TripleSystem


class TripleSystemSecular(TripleSystem):
    id = stateid

    def __init__(self,
                 *args,
                 stabcrit,
                 check_quasi_secular=False,
                 npar=6,
                 **kwargs):

        TripleSystem.__init__(self, *args, **kwargs)

        continue_sim = os.path.isfile(self.outname)
        if not self.mute_output:
            if continue_sim:
                hdfile = h5py.File(self.outname, "r+")
            else:
                hdfile = h5py.File(self.outname, "w")
        else:
            continue_sim = False
            hdfile = None

        self.npar = npar

        self.secularcode = tsunami.Okinami(1, 1)

        # Safety checks
        if (self.gr or self.tide) and self.nocol:
            print("Post-Newtonian/tides enabled but collision disabled, re-enabling collisions")
            self.nocol = False
        if (self.gr or self.tide) and self.R1 == 0.0:
            print("Post-Newtonian/tides enabled but body m1 has zero radius, setting to the Schwartzchild radius")
            self.R1 = 2*self.m1/self.secularcode.speed_of_light**2
        if (self.gr or self.tide) and self.R2 == 0.0:
            print("Post-Newtonian/tides enabled but body m2 has zero radius, setting to the Schwartzchild radius")
            self.R2 = 2*self.m2/self.secularcode.speed_of_light**2

        if continue_sim and not self.restart:
                y, self.k1, self.tau1, self.k2, self.tau2, self.begin_time, hascollided, isdynstable = self.load_secular_run(hdfile)
                y[self.id.g1] = self.KU.mod2pi(y[self.id.g1])
                y[self.id.g2] = self.KU.mod2pi(y[self.id.g2])
                y[self.id.h1] = self.KU.mod2pi(y[self.id.h1])

                self.recordpar = list(self.load_record_par(hdfile))
                print("Restarting from time = {:g} yr".format(self.begin_time))

                self.e1_0 = self.e1 = y[self.id.e1]
                self.e2_0 = self.e2 = y[self.id.e2]
                self.ome1_0 = self.ome1 = y[self.id.g1]
                self.ome2_0 = self.ome2 = y[self.id.g2]
                self.Ome2_0 = self.Ome2 = y[self.id.h1]
                self.a1_0 = self.a1 = y[self.id.a1]
                self.i_mut_0 = self.i_mut = tsunami.Okinami.itot_from_y(y)
                self.a2_0 = self.a2 = y[self.id.a2]
                self.m1, self.m2, self.m3 = y[self.id.m1], y[self.id.m2], y[self.id.m3]
                self.R1, self.R2, self.R3 = y[self.id.R1], y[self.id.R2], y[self.id.R3]
        else:
            if continue_sim and self.restart:
                print("Restarting from initial conditions")
                m1, m2, m3, a1, a2, e1, e2, i_mut, ome1, ome2, Ome2, R1, R2, R3, \
                    k1, k2, k3, tau1, tau2, tau3, nu1, nu2, gyr1, gyr2, gyr3, \
                    spin1, spin2, spin3, obl1, obl2, obl3, azim1, azim2, azim3 = self.load_kepler_ic(hdfile)
                del hdfile[self.hdf_ic_secular]
                del hdfile[self.hdf_time_data]

                self.e1_0 = self.e1 = e1
                self.e2_0 = self.e2 = e2
                self.ome1_0 = self.ome1 = ome1
                self.ome2_0 = self.ome2 = ome2
                self.Ome2_0 = self.Ome2 = Ome2
                self.a1_0 = self.a1 = a1
                self.i_mut_0 = self.i_mut = i_mut
                self.a2_0 = self.a2 = a2
                self.m1, self.m2, self.m3 = m1, m2, m3
                self.R1, self.R2, self.R3 = R1, R2, R3
                self.k1, self.tau1, self.k2, self.tau2 = k1, tau1, k2, tau2
                self.k3, self.tau3 = k3, tau3
                self.gyr1, self.gyr2, self.gyr3 = gyr1, gyr2, gyr3
                self.spin1, self.spin2, self.spin3 = spin1, spin2, spin3
                self.obl1, self.obl2, self.obl3 = obl1, obl2, obl3
                self.azim1, self.azim2, self.azim3 = azim1, azim2, azim3

            y = np.array([self.e1_0, self.e2_0,
                          self.ome1_0, self.ome2_0,
                          self.Ome2_0, self.a1_0, 0.0,
                          self.a2_0, self.m1, self.m2, self.m3,
                          self.R1, self.R2, self.R3])

            tsunami.Okinami.compute_H_in_y(y, self.i_mut_0)
            isdynstable = tsunami.Okinami.is_mardling_stable(self.a1_0, self.m1, self.m2,
                                                             self.a2_0, self.m3, self.e2_0, self.i_mut_0) if stabcrit else True
            hascollided = (self.R1 + self.R2) > (self.a1_0 * (1 - self.e1_0))

            self.begin_time = 0.0
            self.recordpar = []

            if not self.mute_output:
                self.save_secular_run(hdfile,
                                      y, self.k1, self.tau1, self.k2, self.tau2, self.begin_time,
                                      hascollided, isdynstable)
        if not self.mute_output:
            self.save_kepler_ic(hdfile,
                                self.m1, self.m2, self.m3, self.a1, self.a2, self.e1, self.e2, self.i_mut,
                                self.ome1, self.ome2, self.Ome2, self.R1, self.R2, self.R3,
                                self.k1, self.k2, self.k3, self.tau1, self.tau2, self.tau3,
                                None, None, self.gyr1, self.gyr2, self.gyr3, self.spin1, self.spin2, self.spin3,
                                self.obl1, self.obl2, self.obl3, self.azim1, self.azim2, self.azim3)

            hdfile.close()

        self.info_system()

        self.Tzkl = self.secularcode.get_zkl_timescale(self.m1, self.m2, self.a1, self.m3, self.a2,
                                                       self.e2) / self.time_unit
        print("Tzkl = {:e} yr".format(self.Tzkl))

        self.secularcode.initialize_constants(self.m1, self.m2, self.a2_0, self.m3)

        abs_err = 1e-11
        rel_err = 1e-11
        self.secularcode.set_integrator_tolerance(abs_err=abs_err, rel_err=rel_err)
        self.secularcode.check_stability = stabcrit
        self.secularcode.tides = self.tide
        self.secularcode.check_collisions = not self.nocol
        self.secularcode.check_sec = check_quasi_secular

        self.secularcode.gr = self.gr
        self.secularcode.pn1 = not self.nopn1
        self.secularcode.pn2 = not self.nopn2
        self.secularcode.pn25 = not self.nopn25

        self.secularcode.set_tidal_parameters(self.k1, self.tau1, self.k2, self.tau2)

        self.y_lastsane = self.y_0 = self.y = y.copy()

        if np.any(np.isnan(y)):
            print("Found bad values in initial state vector:\ny = {:}\nAborting".format(y))
            self.skip_run = True
        if not isdynstable:
            print("The system is already unstable:\ny = {:}\nAborting".format(y))
            self.skip_run = True
        if hascollided:
            print("The system is already collided:\ny = {:}\nAborting".format(y))
            self.skip_run = True

    def save_secular_run(self, hdfile,
                         y, k1, tau1, k2, tau2, time, hascollided, isdynstable,
                         tstop=None):
        update_hdf5(hdfile, self.hdf_ic_secular + "/y", y)
        update_hdf5(hdfile, self.hdf_ic_secular + "/k1", k1)
        update_hdf5(hdfile, self.hdf_ic_secular + "/tau1", tau1)
        update_hdf5(hdfile, self.hdf_ic_secular + "/k2", k2)
        update_hdf5(hdfile, self.hdf_ic_secular + "/tau2", tau2)
        update_hdf5(hdfile, self.hdf_ic_secular + "/time", time)
        update_hdf5(hdfile, self.hdf_ic_secular + "/hascollided", hascollided, t=bool)
        update_hdf5(hdfile, self.hdf_ic_secular + "/isdynstable", isdynstable, t=bool)
        if tstop is not None: update_hdf5(hdfile, self.hdf_ic_secular + "/tstop", tstop)

    def load_secular_run(self, hdfile):
        y = hdfile[self.hdf_ic_secular + "/y"][...]
        k1 = hdfile[self.hdf_ic_secular + "/k1"][()]
        tau1 = hdfile[self.hdf_ic_secular + "/tau1"][()]
        k2 = hdfile[self.hdf_ic_secular + "/k2"][()]
        tau2 = hdfile[self.hdf_ic_secular + "/tau2"][()]
        time = hdfile[self.hdf_ic_secular + "/time"][...]
        hascollided = hdfile[self.hdf_ic_secular + "/hascollided"][...]
        isdynstable = hdfile[self.hdf_ic_secular + "/isdynstable"][...]
        return y, k1, tau1, k2, tau2, time, hascollided, isdynstable

    def run(self, ftime, dtout, logtime=False, timing=True):
        """
        All times in input here are in yr and are converted into code units afterwards
        :param ftime:
        :param dtout:
        :param logtime:
        :return:
        """
        ftime = (ftime - self.begin_time) * self.time_unit
        dtout = dtout * self.time_unit
        if logtime:
            pass
        if self.skip_run: return

        recordpar = self.recordpar
        self.time = 0.0
        while self.time < ftime:
            self.time = self.time + dtout
            self.secularcode.evolve_system(self.y, self.time)
            self.time = self.secularcode.ctime

            if np.any(np.isnan(self.y)):
                warnings.warn("Found bad values in state vector:\ny = {:}\nSkipping save".format(self.y), Warning)
                break
            else:
                self.y_lastsane = self.y.copy()

            if self.secularcode.hascollided:
                print("\nSystem has collided at time {:g} yr".format(self.time / self.time_unit))
                print("peri {:g} < {:g} R1 + R2".format(self.y[self.id.a1] * (1 - self.y[self.id.e1]), self.y[self.id.R1] + self.y[self.id.R2]))
                break
            if not self.secularcode.isdynstable:
                print("\nSystem is unstable")
                print("a1, e1 = {:g}, {:g}\na2, e2 = {:g}, {:g}".format(self.y[self.id.a1], self.y[self.id.e1], self.y[self.id.a2], self.y[self.id.e2]))
                break
            if self.secularcode.stopsec:
                print("\nSystem entered quasi-secular evolution")
                print("a1, e1 = {:g}, {:g}\na2, e2 = {:g}, {:g}".format(self.y[self.id.a1], self.y[self.id.e1], self.y[self.id.a2], self.y[self.id.e2]))
                break

            if self.killed: break

            self.time = self.time + dtout
            recordpar.append(self.begin_time + self.time / self.time_unit)

            recordpar.append(self.y[self.id.e1])
            recordpar.append(self.y[self.id.e2])
            recordpar.append(self.y[self.id.a1])
            recordpar.append(self.y[self.id.a2])
            recordpar.append(self.secularcode.itot_from_y(self.y))
            if timing: print("Time (secular) = {:05.5g} yr\033[K".format(self.begin_time + self.time / self.time_unit), end="\033[K\r",
                             flush=True)

        print("1-e1max = {:g}, e1min = {:g}".format(1 - self.secularcode.e1max, self.secularcode.e1min))
        self.recordpar = np.array(recordpar)

        if not self.mute_output:
            hdfile = h5py.File(self.outname, "r+")
            self.save_secular_run(hdfile,
                                  self.y_lastsane, self.k1, self.tau1, self.k2, self.tau2,
                                  self.begin_time + self.time / self.time_unit,
                                  self.secularcode.hascollided, self.secularcode.isdynstable,
                                  tstop=self.time)
            self.save_record_par(hdfile, self.recordpar, self.npar)
            hdfile.close()

        print("Peri1, min:", self.a1_0 * (1 - self.secularcode.e1max))
        print("R1+R2:", self.R1 + self.R2)
        print("Final t = {:g} yr".format(self.begin_time + self.time / self.time_unit))


if __name__ == "__main__":
    args = okinami_option_parser().parse_args()

    m1, m2, m3, R1, R2, R3, a1, e1, a2, e2, i_mut, ome1, ome2, Ome2, k1, k2, k3, \
        tau1, tau2, tau3, spin1, spin2, spin3, obl1, obl2, obl3 = process_args(args)
    azim1 = azim2 = azim3 = 0.0

    kzt = TripleSystemSecular(m1, m2, m3,
                              a1, a2,
                              e1, e2,
                              i_mut,
                              ome1, ome2,
                              Ome2,
                              R1, R2, R3,
                              k1, k2, k3,
                              tau1, tau2, tau3,
                              args.gyr1, args.gyr2, args.gyr3,
                              spin1, spin2, spin3, obl1, obl2, obl3, azim1, azim2, azim3,
                              args.nopn1,  args.nopn2, args.nopn25, args.nopn3, args.nopn35,
                              args.nocol, args.tide, args.gr,
                              check_quasi_secular=args.qsec,
                              stabcrit=args.stabcrit,
                              filename=args.filename,
                              mute_output=args.mute_output,
                              restart=args.restart)

    kzt.run(args.ftime, args.dt_out)

    kzt.plot(task=args.task, wSpins=False)
