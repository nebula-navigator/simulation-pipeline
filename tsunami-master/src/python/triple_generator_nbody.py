import tsunami
import numpy as np
import argparse
import h5py
import os
import itertools
from triple_generator_common import tsunami_option_parser, process_args, update_hdf5, TripleSystem


class TripleSystemNbody(TripleSystem):

    def __init__(self,
                 *args,
                 nu1, nu2,
                 begin_time=0.0,
                 invariant=False,
                 npar=12,
                 schw=None,
                 **kwargs,
                 ):

        TripleSystem.__init__(self, *args, **kwargs)

        self.nu1_0 = self.nu1 = nu1
        self.nu2_0 = self.nu2 = nu2
        self.invariant_plane = invariant

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

        self.Rschw_mult = schw
        if self.Rschw_mult is not None:
            print("Setting all radii {:g} times the Schwartzchild radius".format(self.Rschw_mult))
            self.R1 = self.Rschw_mult * 2 * self.m1 / self.KU.speed_of_light ** 2
            self.R2 = self.Rschw_mult * 2 * self.m2 / self.KU.speed_of_light ** 2
            self.R3 = self.Rschw_mult * 2 * self.m3 / self.KU.speed_of_light ** 2

        # Safety checks
        if (self.gr or self.tide) and self.nocol:
            print("Post-Newtonian/tides enabled but collision disabled, re-enabling collisions")
            self.nocol = False
        if (self.gr or self.tide) and self.R1 == 0.0:
            print("Post-Newtonian/tides enabled but body m1 has zero radius, setting to the Schwartzchild radius")
            self.R1 = 2 * self.m1 / self.KU.speed_of_light ** 2
        if (self.gr or self.tide) and self.R2 == 0.0:
            print("Post-Newtonian/tides enabled but body m2 has zero radius, setting to the Schwartzchild radius")
            self.R2 = 2 * self.m2 / self.KU.speed_of_light ** 2
        if (self.gr or self.tide) and self.R3 == 0.0:
            print("Post-Newtonian/tides enabled but body m3 has zero radius, setting to the Schwartzchild radius")
            self.R3 = 2 * self.m3 / self.KU.speed_of_light ** 2

        # If even terms are actually enabled, useful for PN-corrected elements
        self.pn1 = self.gr and not self.nopn1
        self.pn2 = self.gr and not self.nopn2
        self.pn3 = self.gr and not self.nopn3

        if continue_sim and not self.restart:
            m, p, v, s, kaps, taus, rads, gyrs, self.begin_time, collision, breakup = self.load_nbody_run(hdfile)
            m1, m2, m3, a1, a2, e1, e2, i_mut, ome1, ome2, Ome2, R1, R2, R3, \
                k1, k2, k3, tau1, tau2, tau3, nu1, nu2, gyr1, gyr2, gyr3, \
                spin1, spin2, spin3, obl1, obl2, obl3, azim1, azim2, azim3 = self.load_kepler_ic(hdfile)
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

            self.recordpar = list(self.load_record_par(hdfile))
            self.m1 = m[0]
            self.m2 = m[1]
            self.m3 = m[2]
            self.R1, self.R2, self.R3 = rads[0], rads[1], rads[2]
            self.k1, self.k2, self.k3 = kaps[0], kaps[1], kaps[2]
            self.tau1, self.tau2, self.tau3 = taus[0], taus[1], taus[2]
            self.gyr1, self.gyr2, self.gyr3 = gyrs[0], gyrs[1], gyrs[2]
            print("Restarting from time = {:g} yr".format(self.begin_time))
        else:
            if continue_sim and self.restart:
                print("Restarting from initial conditions")
                m1, m2, m3, a1, a2, e1, e2, i_mut, ome1, ome2, Ome2, R1, R2, R3, \
                    k1, k2, k3, tau1, tau2, tau3, nu1, nu2, gyr1, gyr2, gyr3, \
                    spin1, spin2, spin3, obl1, obl2, obl3, azim1, azim2, azim3 = self.load_kepler_ic(hdfile)
                del hdfile[self.hdf_ic_nbody]
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

            m, p, v = self.create_nbody_ic(self.m1, self.m2, self.m3,
                                           self.a1, self.a2,
                                           self.e1, self.e2,
                                           self.i_mut,
                                           self.ome1, self.ome2,
                                           self.Ome2,
                                           self.nu1, self.nu2,
                                           invariant_plane=self.invariant_plane)
            s = self.setup_spin_vectors(self.spin1, self.spin2, self.spin3, self.obl1, self.obl2, self.obl3,
                                        self.azim1, self.azim2, self.azim3)

            self.begin_time = begin_time
            self.recordpar = []
            kaps = np.array([self.k1, self.k2, self.k3])
            taus = np.array([self.tau1, self.tau2, self.tau3])
            rads = np.array([self.R1, self.R2, self.R3])
            gyrs = np.array([self.gyr1, self.gyr2, self.gyr3])

            if not self.mute_output:
                self.save_nbody_run(hdfile,
                                    m, p, v, s, kaps, taus, rads, gyrs, self.begin_time)

        if not self.mute_output:
            self.save_kepler_ic(hdfile,
                                self.m1, self.m2, self.m3, self.a1, self.a2, self.e1, self.e2, self.i_mut,
                                self.ome1, self.ome2, self.Ome2, self.R1, self.R2, self.R3,
                                self.k1, self.k2, self.k3, self.tau1, self.tau2, self.tau3,
                                self.nu1, self.nu2, self.gyr1, self.gyr2, self.gyr3, self.spin1, self.spin2, self.spin3,
                                self.obl1, self.obl2, self.obl3, self.azim1, self.azim2, self.azim3)

            hdfile.close()

        self.info_system()
        print("ICs in the invariant plane:", self.invariant_plane)

        st = np.array([1, 2, 3])

        self.directcode = tsunami.Tsunami()

        self.directcode.Conf.wEqTides = self.tide
        self.directcode.Conf.wPNs = self.gr
        self.directcode.Conf.pn1 = not self.nopn1
        self.directcode.Conf.pn2 = not self.nopn2
        self.directcode.Conf.pn25 = not self.nopn25
        self.directcode.Conf.pn3 = not self.nopn3
        self.directcode.Conf.pn35 = not self.nopn35

        if self.nocol: self.directcode.Conf.dcoll = 0.0

        self.directcode.add_particle_set(p, v, m, rads, st, s)

        self.directcode.initialize_tidal_parameters(kaps, taus, [1.5, 1.5, 1.5], gyrs)

        self.p, self.v, self.m, self.s = p, v, m, s
        self.kaps, self.taus, self.rads, self.gyrs = kaps, taus, rads, gyrs
        self.minn = self.m1 + self.m2
        self.mout = self.m1 + self.m2 + self.m3

        self.directcode.sync_internal_state(p, v)

    def create_nbody_ic(self, m1, m2, m3,
                        a1, a2,
                        e1, e2,
                        i_mut,
                        ome1, ome2,
                        Ome2,
                        nu1, nu2, invariant_plane=False):

        i1, i2 = tsunami.Okinami.i1_i2_from_itot(m1, m2, m3,
                                                 a1, a2,
                                                 e1, e2,
                                                 i_mut)
        Ome1 = Ome2 + np.pi
        Ome1 = self.KU.mod2pi(Ome1)

        if not invariant_plane:
            i1 = Ome1 = 0.0
            i2 = i_mut

        pos_vel2 = np.array([0., 0., 0., 0., 0., 0.])
        minn = m1 + m2
        pos_vel1 = self.KU.kepl_to_cart(pos_vel2, m1, m2, a1, e1, i1, ome1, Ome1, nu1, pn1=self.pn1, pn2=self.pn2)

        pos_vel_inn_com = (pos_vel1 * m1 + pos_vel2 * m2) / minn
        pos_vel1 -= pos_vel_inn_com
        pos_vel2 -= pos_vel_inn_com

        mout = m1 + m2 + m3
        pos_vel3 = self.KU.kepl_to_cart(np.zeros(6, dtype=np.float64), minn, m3, a2, e2, i2, ome2, Ome2, nu2, pn1=self.pn1, pn2=self.pn2)

        triple_com = (pos_vel1 * m1 + pos_vel2 * m2 + pos_vel3 * m3) / mout
        pos_vel1 -= triple_com
        pos_vel2 -= triple_com
        pos_vel3 -= triple_com

        m = np.array([m1, m2, m3], dtype=np.float64)
        p = np.array([pos_vel1[:3], pos_vel2[:3], pos_vel3[:3]], dtype=np.float64)
        v = np.array([pos_vel1[3:], pos_vel2[3:], pos_vel3[3:]], dtype=np.float64)

        return m, p, v

    def save_nbody_run(self, hdfile,
                       m, p, v, s, kaps, taus, rads, gyrs, time, collision=False, breakup=False):
        update_hdf5(hdfile, self.hdf_ic_nbody + "/m", m)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/p", p)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/v", v)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/s", s)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/kaps", kaps)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/taus", taus)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/rads", rads)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/gyrs", gyrs)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/time", time)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/collision", collision, t=bool)
        update_hdf5(hdfile, self.hdf_ic_nbody + "/breakup", breakup, t=bool)

    def load_nbody_run(self, hdfile):
        m = hdfile[self.hdf_ic_nbody + "/m"][...]
        p = hdfile[self.hdf_ic_nbody + "/p"][...]
        v = hdfile[self.hdf_ic_nbody + "/v"][...]
        s = hdfile[self.hdf_ic_nbody + "/s"][...]
        kaps = hdfile[self.hdf_ic_nbody + "/kaps"][...]
        taus = hdfile[self.hdf_ic_nbody + "/taus"][...]
        rads = hdfile[self.hdf_ic_nbody + "/rads"][...]
        gyrs = hdfile[self.hdf_ic_nbody + "/gyrs"][...]
        time = hdfile[self.hdf_ic_nbody + "/time"][...]
        collision = hdfile[self.hdf_ic_nbody + "/collision"][...]
        breakup = hdfile[self.hdf_ic_nbody + "/breakup"][...]
        return m, p, v, s, kaps, taus, rads, gyrs, time, collision, breakup

    def run(self, ftime, dtout, timing=True):
        print("Final time: {:g} yr  ({:g} P1, {:g} P2)".format(ftime, ftime/self.P1*self.time_unit, ftime/self.P2*self.time_unit))
        print("Output steps: {:g} yr  ({:g} P1, {:g} P2)".format(dtout, dtout/self.P1*self.time_unit, dtout/self.P2*self.time_unit))

        if self.directcode.check_collision():
            id1, id2 = self.directcode.get_collision_indices()
            dp = self.p[id1]-self.p[id2]
            r = np.sqrt((dp*dp).sum())
            R = self.rads[id1] + self.rads[id2]
            print("System is already in collision")
            print("\tParticles ({:d},{:d}) have distance {:g} au\n\tBut the sum of their radii is {:g}".format(id1, id2, r*self.KU.Lscale, R*self.KU.Lscale))
            print("\tSkipping evolution")
            return

        ftime = (ftime - self.begin_time) * self.time_unit
        dtout = dtout * self.time_unit
        recordpar = list(self.recordpar)
        self.time = self.directcode.time = 0.0  #FIXME
        breakup = collision = False
        while self.time < ftime:
            self.time = self.time + dtout
            self.directcode.evolve_system(self.time)
            self.time = self.directcode.time
            self.directcode.sync_internal_state(self.p, self.v, self.s)

            self.orbpars_pvs(self.p, self.v, self.s)
            if self.killed: break

            recordpar.append(self.begin_time + self.time / self.time_unit)
            recordpar.append(self.orb1[1])  # e1
            recordpar.append(self.orb2[1])  # e2
            recordpar.append(self.orb1[0])  # a1
            recordpar.append(self.orb2[0])  # a2
            recordpar.append(self.i_mut)    # i_mut (rad)
            recordpar.append(self.orb1[3])  # ome1 (rad)
            recordpar.append(self.orb2[3])  # ome2 (rad)
            recordpar.append(self.spin1 * self.time_unit)  # spin1 (rad/yr)
            recordpar.append(self.spin2 * self.time_unit)  # spin2 (rad/yr)
            recordpar.append(self.obl1)  # obl1 (rad)
            recordpar.append(self.obl2)  # obl2 (rad)

            if self.directcode.check_collision() > 0:
                print("\nCollision detected at time {:05.5g} yr".format(self.begin_time + self.time / self.time_unit))
                collision = True
                id1, id2 = self.directcode.get_collision_indices()
                break

            breakup = self.check_stability()
            if breakup:
                print("\nBreakup detected at time {:05.5g} yr".format(self.begin_time + self.time / self.time_unit))
                break

            if timing: print("Time (direct) = {:05.5g} yr\033[K".format(self.begin_time + self.time / self.time_unit), end="\033[K\r")

        self.recordpar = np.array(recordpar)
        self.begin_time = self.time

        if not self.mute_output:
            hdfile = h5py.File(self.outname, "r+")
            self.save_record_par(hdfile,
                                 self.recordpar, self.npar)
            self.save_nbody_run(hdfile,
                                self.m, self.p, self.v, self.s,
                                self.kaps, self.taus, self.rads, self.gyrs, self.begin_time + self.time / self.time_unit,
                                collision, breakup)
            hdfile.close()

    def orbpars_pvs(self, p, v, s):
        p_inn = p[0]-p[1]
        v_inn = v[0]-v[1]
        self.orb1 = self.KU.cart_to_kepl(p_inn, v_inn, self.m[0], self.m[1], pn1=self.pn1, pn2=self.pn2, pn3=self.pn3, memmcorr=self.gr)
        p_com_inn = (p[0] * self.m1 + p[1] * self.m2) / self.minn
        v_com_inn = (v[0] * self.m1 + v[1] * self.m2) / self.minn

        p_out = p[2] - p_com_inn
        v_out = v[2] - v_com_inn
        self.orb2 = self.KU.cart_to_kepl(p_out, v_out, self.minn, self.m3, pn1=self.pn1, pn2=self.pn2, pn3=self.pn3, memmcorr=self.gr)

        cosi = np.cos(self.orb1[2])*np.cos(self.orb2[2]) + np.sin(self.orb1[2])*np.sin(self.orb2[2]) * np.cos(self.orb1[4] - self.orb2[4])

        h_vec = np.cross(p_inn, v_inn)
        h = (h_vec*h_vec).sum()**0.5
        h_hat = h_vec/h


        if self.directcode.Conf.wSpins:
            self.spin1 = (s[0]*s[0]).sum()**0.5
            self.spin2 = (s[1]*s[1]).sum()**0.5
            spin1_hat = s[0]/self.spin1
            spin2_hat = s[1]/self.spin2
            self.obl1 = np.arccos((spin1_hat * h_hat).sum())
            self.obl2 = np.arccos((spin2_hat * h_hat).sum())

        self.i_mut = np.arccos(cosi)

    def orbpars_pv2(self, p, v):
        pv = np.hstack((p, v))

        self.orblist = []
        for (i, j), k in zip(itertools.combinations(range(3), r=2), range(2, -1, -1)):
            pv_inn = pv[i] - pv[j]
            minn = self.m[i] + self.m[j]
            orb1 = self.KU.cart_to_kepl(pv_inn, minn)
            com_inn = (pv[i] * self.m[i] + pv[j] * self.m[j]) / minn

            pv_out = pv[k] - com_inn
            mout = minn + self.m[k]
            orb2 = self.KU.cart_to_kepl(pv_out, mout)
            cosi = np.cos(orb1[2])*np.cos(orb2[2]) + np.sin(orb1[2])*np.sin(orb2[2]) * np.cos(orb1[4] - orb2[4])

            i_mut = np.arccos(cosi)
            i_mut = np.degrees(i_mut)
            ene1 = - self.m[i] * self.m[j] / (2*orb1[0])

            self.orblist.append((ene1, orb1,orb2,(i,j,k)))
            if (i,j) == (0,1):
                self.i_mut = i_mut
                self.orb1 = orb1
                self.orb2 = orb2

    def check_stability(self):
        if (self.orb1[0] < 0) or (self.orb2[0] < 0):
            return True
        else:
            return False

    def check_stability2(self):
        if self.orb1[0] > 0 and self.orb2[0] > 0:
            # Still in original stable config
            return False
        self.orblist.sort(key=lambda x: x[0])
        mostb = self.orblist[0]
        ene1, orb1, orb2, (i,j,k) = mostb
        if ene1 > 0:
            # All are unbound
            # Still doing dynamical interactions
            return False
        if ene1 < 0 and orb2[0] < 0:
            # Inner bound, outer unbound
            p_inn = (self.p[i] * self.m[i] + self.p[j] * self.m[j]) / (self.m[i] + self.m[j])
            relp = p_inn - self.p[k]
            d = (relp*relp).sum()**0.5
            if d > 10 * orb1[0]:
                # Binary single is far, stop here
                return True
        return False

        #if (self.orb1[0] < 0) or (self.orb2[0] < 0):
        #    return True
        #else: return False

    def save_final_info(self, hdfile, breakup, collision):
        lastorb = self.orblist[0]
        abin = lastorb[1][0]
        ebin = lastorb[1][1]
        i, j = lastorb[3][0], lastorb[3][1]
        m1, m2 = self.m[i], self.m[j]

        update_hdf5(hdfile, "data_nb/breakup", breakup, t=bool)
        update_hdf5(hdfile, "data_nb/collision", collision, t=bool)
        update_hdf5(hdfile, "data_nb/abin", abin)
        update_hdf5(hdfile, "data_nb/ebin", ebin)
        update_hdf5(hdfile, "data_nb/i", i, t=np.int)
        update_hdf5(hdfile, "data_nb/j", j, t=np.int)
        update_hdf5(hdfile, "data_nb/m1", m1)
        update_hdf5(hdfile, "data_nb/m2", m2)


if __name__ == "__main__":
    args = tsunami_option_parser().parse_args()

    m1, m2, m3, R1, R2, R3, a1, e1, a2, e2, i_mut, ome1, ome2, Ome2, k1, k2, k3, \
        tau1, tau2, tau3, spin1, spin2, spin3, obl1, obl2, obl3 = process_args(args)
    nu1, nu2 = np.radians([args.nu1, args.nu2])
    azim1 = azim2 = azim3 = 0.0
    kzt = TripleSystemNbody(m1, m2, m3,
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
                            args.nopn1, args.nopn2, args.nopn25, args.nopn3, args.nopn35,
                            args.nocol, args.tide, args.gr,
                            nu1=nu1, nu2=nu2,
                            filename=args.filename,
                            mute_output=args.mute_output,
                            invariant=args.invariant,
                            restart=args.restart)
    kzt.run(args.ftime, args.dt_out)
    kzt.plot(task=args.task, wSpins=kzt.directcode.Conf.wSpins)
