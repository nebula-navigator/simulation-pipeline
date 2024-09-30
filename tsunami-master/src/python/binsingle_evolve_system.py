"""
binsingle_evolve_system.py

This script runs an individual binary-single interaction, analyzes its evolution and returns the outcome parameters.

It can be run as a standalone script, or used in together with binsingle_initial_conditions.py and binsingle_run.py
In the latter case, it will be run automatically by binsingle_run.py, when called by the <setname>.sh bash script

The outcome parameters of the simulation run are encapsulated in the StatusIsol class. An instance of this class is
returned at the end of a simulation, from the function run_system
Its data members are:
    * begin_time :  time (in N-body units) when the binary-single begin interacting
    * end_time : time (in N-body units) when the binary-single stop interacting
    * ibin : tuple with the particle indices of the (final) escaping binary
    * isin : index of the (final) single
    * Nexcurs : number of excursions
    * Nexchan : number of exchanges
    * collision : boolean, was there a collision?
    * collid : indices of collided particles
    * orb_inn : keplerian orbital elements of the final (bound) elliptical binary
    * orb_out : keplerian orbital elements of the outer (unbound) hyperbolic binary
    * ex_tbegin : begin time of the last excursion
    * status : string containing the status of the interaction. Can be:
        - "begin" : beginning of the simulation, before any interaction
        - "excursion" : the triple is undergoing a bound excursion
        - "breakup" : the triple has broken up into 3 unbound particles
        - "interacting" : the triple is still interacting
        - "escape" : the triple has decayed into an unbound binary-single on diverging orbits
        - "will_interact" : the triple is in an unbound binary-single, but it will interact in the future (i.e. converging orbits)

The simulation is run until either the final time (ftime) is reached, or when the triple has reached its final state
(an unbound, diverging binary-single). The status of the simulation is tracked by the stopping_condition function, which
will stop the simulation as soon as the triple decays into a binary and a single (or three unbound bodies, untested)
"""

import tsunami
import numpy as np
import argparse
import os
import h5py
import pandas as pd
from collections import deque

KU = tsunami.KeplerUtils()

ib1 = 0
ib2 = 1
isi = 2
ind_loop = {0: [0, 1, 2], 1: [0, 2, 1], 2: [1, 2, 0]}
ind_arr = np.array([0, 1, 2])
coll_folder_suff = "_coll"


def setup_system(m1, m2, m3, R1, R2, R3, a_inn, e_inn, ome_inn, nu_inn, vinf, b, i_out, ome_out, Ome_out,
                 d_sin_a_inn=150):
    pos_vel2 = np.array([0., 0., 0., 0., 0., 0.])
    pos_vel1 = KU.kepl_to_cart(pos_vel2, m1, m2, a_inn, e_inn, 0.0, ome_inn, 0.0, nu_inn)

    minn = m1 + m2
    pos_vel_inn_com = (pos_vel1 * m1 + pos_vel2 * m2) / minn

    d_sin = d_sin_a_inn * a_inn

    mout = minn + m3
    a_sin = - mout / vinf ** 2

    e_sin = (1. + (b / a_sin) ** 2) ** 0.5
    if e_sin == 1:
        # e_sin = (1 + 0.5*(b / a_sin)**2)
        latum = - b * b / a_sin
        e_sin = 1 + 1e-9
        rp = latum / (1 + e_sin)
        a_sin = rp / (1 - e_sin)
    else:
        latum = a_sin * (1. - e_sin ** 2)
        rp = a_sin * (1. - e_sin)
    if d_sin < rp:
        raise ValueError("Initial distance is smaller than pericenter distance, {:e} < {:e}".format(d_sin, rp))

    nu_sin = latum / d_sin - 1
    nu_sin = nu_sin / e_sin
    nu_sin = - np.arccos(nu_sin)

    E_sin = KU.nu_to_E(nu_sin, e_sin)
    M_sin = KU.E_to_M(E_sin, e_sin)

    Pout = 2 * np.pi * (np.fabs(a_sin) ** 3 / mout) ** 0.5
    time_to_enc = np.fabs(M_sin) / (2 * np.pi) * Pout
    Pinn = 2 * np.pi * (a_inn ** 3 / minn) ** 0.5

    pos_vel3 = KU.kepl_to_cart(pos_vel_inn_com, minn, m3, a_sin, e_sin, i_out, ome_out, Ome_out, nu_sin)

    triple_com = (pos_vel1 * m1 + pos_vel2 * m2 + pos_vel3 * m3) / mout
    pos_vel1 -= triple_com
    pos_vel2 -= triple_com
    pos_vel3 -= triple_com

    pos_vel_list = np.array([pos_vel1, pos_vel2, pos_vel3])

    # Add here GR corrections
    m = np.array([m1, m2, m3], dtype=np.float64)
    R = np.array([R1, R2, R3], dtype=np.float64)
    p = pos_vel_list[:, 0:3].copy(order='C')  # Ensures that vectors are contiguous
    v = pos_vel_list[:, 3:6].copy(order='C')

    dp_inn = p[ib1] - p[ib2]
    dv_inn = v[ib1] - v[ib2]
    orb_inn = KU.cart_to_kepl(dp_inn, dv_inn, m1, m2)
    pcom_inn = p[ib1] * m[ib1] + p[ib2] * m[ib2]
    pcom_inn /= minn
    vcom_inn = v[ib1] * m[ib1] + v[ib2] * m[ib2]
    vcom_inn /= minn
    orb_out = KU.cart_to_kepl(p[isi] - pcom_inn, v[isi] - vcom_inn, minn, m3)

    print("m1, m2, m3 = {:3.2f}, {:3.2f}, {:3.2f}".format(m1, m2, m3))
    print("P_inn [yr/2pi] = {:3.2f}".format(Pinn))
    print("a_inn [au], e_inn, i_inn = {:3.2f}, {:1.2g}, {:1.2g}".format(orb_inn[0], orb_inn[1], np.degrees(orb_inn[2])))
    print("ome_inn, Ome_inn, nu_inn = {:3.2f}, {:3.2f}, {:3.2f}".format(np.degrees(orb_inn[3]), np.degrees(orb_inn[4]),
                                                                        np.degrees(orb_inn[5])))
    print("P_out [yr/2pi] = {:3.2f}".format(Pout))
    print("a_out [au], e_out, i_out = {:3.2f}, {:1.2g}, {:1.2g}".format(orb_out[0], orb_out[1], np.degrees(orb_out[2])))
    print("ome_out, Ome_out, nu_out = {:3.2f}, {:3.2f}, {:3.2f}".format(np.degrees(orb_out[3]), np.degrees(orb_out[4]),
                                                                        np.degrees(orb_out[5])))
    v_to_kms = 2 * np.pi * KU.au2km / KU.yr2sec
    print("Vinf [km/s] = {:3.2f}".format(vinf * v_to_kms))
    print("B [au], B/a_bin = {:3.2f}, {:1.2g}".format(b, b / a_inn))

    print("time_to_enc [yr/2pi] = {:3.2f}".format(time_to_enc))
    print("time_to_enc/P_inn = {:1.2g}".format(time_to_enc / Pinn))

    return p, v, m, R, orb_inn, orb_out, time_to_enc


def run_system(p, v, m, R, ftime, dtout, time_to_enc, save_traj=False, logging=False, timing=False,
               fname="testrun", setname=None, gw=True, collisions=True, exmap=False):
    code = tsunami.Tsunami()
    code.Conf.wPNs = gw
    if not collisions and not gw:
        R[:] = 0.0
    elif not collisions and gw:
        print("WARNING: cannot disable collisions without disabling PN terms\nadd the -nogw options")

    coll_folder = setname + coll_folder_suff if setname is not None else ""
    logname = setname + "_" + fname + "_log.txt" if setname is not None else fname + "_log.txt"

    buffer_length = 100
    time_buffer = deque(maxlen=buffer_length)
    coord_buffer = deque(maxlen=buffer_length)

    code.add_particle_set(p, v, m, R, np.zeros_like(m, dtype=np.int64))
    code.sync_internal_state(p, v)

    logfile = None if not logging else logname
    status = StatusIsol(logfile=logfile, exmap=exmap)

    time = 0.0
    totp, totv, tott = [p.copy()], [v.copy()], [time]
    while True:
        time = time + dtout if not exmap else time * (1 + 1e-15)

        code.evolve_system(time)
        time = code.time
        code.sync_internal_state(p, v)

        time_buffer.append(time)
        coord_buffer.append((p, v))

        if save_traj:
            totp.append(p.copy())
            totv.append(v.copy())
            tott.append(time)

        stop = stopping_condition(m, p, v, time, status)
        if stop and not code.check_collision() > 0 and not time < time_to_enc:
            break

        # Check for collisions, stop evolution
        if code.check_collision() > 0:
            print("\nCollision at time: {:g} yr".format(time))
            # index of colliding particles, always in ascending order
            id1, id2 = code.get_collision_indices()

            status.collision = True
            status.collid = [id1, id2]
            status.end_time = time
            nk = np.delete(ind_arr, [id1, id2])[0]
            status.orb_inn = KU.cart_to_kepl(p[id1] - p[id2], v[id1] - v[id2], m[id1], m[id2])
            mc = m[id1] + m[id2]
            pc_com, vc_com = p[id1] * m[id1] + p[id2] * m[id2], v[id1] * m[id1] + v[id2] * m[id2]
            pc_com, vc_com = pc_com / mc, vc_com / mc
            status.orb_out = KU.cart_to_kepl(pc_com - p[nk], vc_com - v[nk], mc, m[nk])

            # code.save_restart_file(os.path.join(coll_folder, fname + ".bin"))
            break

        if time > ftime:
            # Timeout
            # code.save_restart_file(os.path.join(coll_folder, fname+"_restart.bin"))
            if status.status == "excursion":
                i, j = status.ibin
                k = status.isin
                status.orb_inn = KU.cart_to_kepl(p[i] - p[j], v[i] - v[j], m[i], m[j])
                mt = m[i] + m[j]
                pc_com, vc_com = p[i] * m[i] + p[j] * m[j], v[i] * m[i] + v[j] * m[j]
                pc_com, vc_com = pc_com / mt, vc_com / mt
                status.orb_out = KU.cart_to_kepl(pc_com - p[k], vc_com - v[k], mt, m[k])
            else:
                status.orb_inn = status.orb_out = np.repeat(-1, 6)
            break

        if timing: print("time={:2.2e}/{:2.2e}, d_sin/a_bin={:e}".format(time, ftime, status.d_sin_inn_ratio),
                         end="\033[K\r")

    if save_traj:
        totp = np.vstack(totp)
        totv = np.vstack(totv)
        tott = np.array(tott)
        h5f = h5py.File(fname + "_traj.h5", "w")
        # All n-body units
        h5f["p"] = totp
        h5f["v"] = totv
        h5f["t"] = tott
        h5f.close()

    if status.collision:
        h5f = h5py.File(os.path.join(coll_folder, fname + "_snapshots.h5"), "w")
        h5f["id1"] = status.collid[0]
        h5f["id2"] = status.collid[1]
        h5f["mass"] = m
        h5f["rad"] = R
        for it in range(len(time_buffer)):
            lastit = len(time_buffer) - it - 1
            g = h5f.create_group("{:03d}".format(it))
            g["time"] = time_buffer[lastit]
            g["pos"] = coord_buffer[lastit][0]
            g["vel"] = coord_buffer[lastit][1]
        h5f.close()

    return status


class StatusIsol:
    def __init__(self, logfile=None, exmap=False):
        self.status = "begin"  # status of the interaction
        self.begin_time = -1.0  # when the binary-single begin interacting
        self.end_time = -1.0  # when the binary-single stop interacting
        self.ibin = (-1, -1)  # indices of the (final) binary
        self.isin = -1  # indices of the (final) single
        self.Nexcurs = 0  # number of excursions
        self.Nexchan = 0  # number of exchanges
        self.collision = False  # was there a collision?
        self.collid = [-1, -1]  # indices of collided particles
        self.orb_inn = None  # keplerian elements of the final (bound) elliptical binary
        self.orb_out = None  # keplerian elements of the outer (unbound) hyperbolic binary
        self.ex_tbegin = -1  # begin time of the last excursion

        self.ex_vdotr = None
        self.ex_list = []
        self.ex_last = None

        self.d_sin_inn_ratio = -1.0
        self.ex_Pinn_Nmean = (-1, -1)

        self.logfile = logfile if logfile is None else open(logfile, 'w')
        self.exmap = exmap

    def assign_binary_and_check_exchange(self, i, j, k):
        exchange = self.ibin != (i, j) if self.ibin != (-1, -1) else False
        self.ibin = i, j
        self.isin = k
        if exchange: self.Nexchan += 1

    def average_Pinn(self, Pinn):
        avePinn, Nmean = self.ex_Pinn_Nmean
        if Nmean == -1:
            self.ex_Pinn_Nmean = Pinn, 1
        else:
            new_avePinn = (avePinn * Nmean + Pinn) / (Nmean + 1)
            self.ex_Pinn_Nmean = new_avePinn, Nmean + 1

    def __repr__(self):
        repstr = "Status = {:s}\nBegin time = {:g} yr/2pi\n".format(self.status, self.begin_time)
        repstr = repstr + "End time = {:g} yr/2pi\nBinary = {:d},{:d}\n".format(self.end_time, *self.ibin)
        repstr = repstr + "Single = {:d}\nNexcurs = {:d}\nNexchan = {:d}\n".format(self.isin, self.Nexcurs,
                                                                                   self.Nexchan)
        repstr = repstr + "Collision = {:}\nColl ids = {:d},{:d}\n".format(self.collision, *self.collid)
        repstr = repstr + "Final orbits, a [au], e, i [rad], ome [rad], Ome [rad], nu [rad]\n"
        repstr = repstr + "Inner = {:g}, {:g}, {:g}, {:g}, {:g}, {:g}\n".format(*self.orb_inn)
        repstr = repstr + "Outer = {:g}, {:g}, {:g}, {:g}, {:g}, {:g}".format(*self.orb_out)
        return repstr


class Excursion:
    def __init__(self, mbin, msin, ainn, einn, iinn, Omeinn,
                 aout, eout, iout, Omeout, tapo):
        self.mbin = mbin
        self.msin = msin
        self.ainn = ainn
        self.einn = einn
        self.iinn = iinn
        self.Omeinn = Omeinn
        self.aout = aout
        self.eout = eout
        self.iout = iout
        self.Omeout = Omeout
        self.tapo = tapo

    def __repr__(self):
        repstr = "tapo = {:g} yr/2pi, mbin = {:g}, msin = {:g}\n".format(self.tapo, self.mbin, self.msin)
        repstr = repstr + "ainn = {:g} au, einn = {:g}, iinn = {:g}, Oinn = {:g}\n".format(self.ainn, self.einn,
                                                                                           np.degrees(self.iinn),
                                                                                           np.degrees(self.Omeinn))
        repstr = repstr + "aout = {:g} au, eout = {:g}, iout = {:g}, Oout = {:g}".format(self.aout, self.eout,
                                                                                         np.degrees(self.iout),
                                                                                         np.degrees(self.Omeout))
        return repstr

    def as_dict(self):
        return {'mbin': self.mbin, 'msin': self.msin, 'ainn': self.ainn, 'einn': self.einn, 'iinn': self.iinn,
                'Omeinn': self.Omeinn,
                'aout': self.aout, 'eout': self.eout, 'iout': self.iout, 'Omeout': self.Omeout, 'tapo': self.tapo}

    def as_string(self, nex=0, fname=""):
        print(fname, nex)
        string = "{:s} {:d} {:1.16e} {:1.16e} {:1.16e}".format(fname, nex, self.tapo, self.mbin, self.msin)
        string = string + " {:1.16e} {:1.16e} {:1.16e} {:1.16e}".format(self.ainn, self.einn, self.iinn, self.Omeinn)
        string = string + " {:1.16e} {:1.16e} {:1.16e} {:1.16e}".format(self.aout, self.eout, self.iout, self.Omeout)
        return string


def stopping_condition(m, p, v, time, sta: StatusIsol, a_escape_mult=20, time_escape_mult=20):
    inva_iter = np.zeros(3)
    mutual_dist = np.zeros(3)
    for key, values in ind_loop.items():
        ip, jp = values[:2]
        mtot_iter = m[ip] + m[jp]
        dp_iter = p[ip] - p[jp]
        dv_iter = v[ip] - v[jp]
        inva_iter[key], mutual_dist[key] = inv_semimajor_d(dp_iter, dv_iter, mtot_iter)

    imax = inva_iter.argmax()
    inva_max = inva_iter[imax]
    i, j, k = ind_loop[imax]
    inva_ik, inva_jk = inva_iter[np.delete(ind_arr, imax)]
    was_excursion = sta.status == "excursion"

    if inva_max < 0:
        # Not bound
        if inva_ik < 0 and inva_ik < 0:
            # Nothing is bound, let's check if bound to others' com
            inva_2com = np.zeros(3)
            for key, values in ind_loop.items():
                ip, jp, kp = values
                p2com = p[ip] * m[ip] + p[jp] * m[jp]
                v2com = v[ip] * m[ip] + v[jp] * m[jp]
                m2com = m[ip] + m[jp]
                p2com /= m2com
                v2com /= m2com

                dp_2com = p2com - p[kp]
                dv_2com = v2com - v[kp]
                inva_2com[key], _ = inv_semimajor_d(dp_2com, dv_2com, m2com + m[k])
            bound2com = inva_2com > 0
            if not bound2com.all():
                # Really nothing is bound, let's check the mutual distances
                if sta.status != "breakup":
                    if sta.logfile: sta.logfile.write("t={:1.3e} : ({:d},{:d},{:d}) : no bound particles, "
                                                      "breakup candidate\n".format(time, i, j, k))
                    sta.end_time = time
                    sta.status = "breakup"
                else:
                    Etot = total_energy(p, v, m)
                    mtot = m.sum()
                    t_dyn = 0.5 * (mtot ** 5 / (np.fabs(Etot) ** 3)) ** 0.5
                    t_elapse = time - sta.end_time
                    time_has_elapsed = t_elapse > t_dyn * time_escape_mult
                    if time_has_elapsed:
                        # It's a breakup
                        if sta.logfile: sta.logfile.write("t={:1.3e} : breakup completed\n".format(time))
                        orbs2c = []
                        for key, values in ind_loop.items():
                            ip, jp, kp = values
                            p2com = p[ip] * m[ip] + p[jp] * m[jp]
                            v2com = v[ip] * m[ip] + v[jp] * m[jp]
                            m2com = m[ip] + m[jp]
                            p2com /= m2com
                            v2com /= m2com

                            dp_2com = p2com - p[kp]
                            dv_2com = v2com - v[kp]
                            orbs2c.append(KU.cart_to_kepl(dp_2com, dv_2com, m2com, m[kp]))

                        sta.orb_inn = orbs2c[0]
                        sta.orb_out = orbs2c[1]
                        return True
            else:
                # They are still interacting
                if sta.status != "interacting":
                    if sta.logfile: sta.logfile.write("t={:1.3e} : [{:d},{:d},{:d}] : no bound particles, "
                                                      "but particles are interacting\n".format(time, i, j, k))
                    sta.status = "interacting"
        else:
            # They are still interacting
            if sta.status != "interacting":
                if sta.logfile: sta.logfile.write("t={:1.3e} : [{:d},{:d},{:d}] : no bound particles, "
                                                  "but particles are interacting\n".format(time, i, j, k))
                sta.status = "interacting"

    if inva_max > 0:
        # i, j are bound binary, k is third body
        # ibin, isin are set only once per classification
        pcom_inn = p[j] * m[j] + p[i] * m[i]
        vcom_inn = v[j] * m[j] + v[i] * m[i]
        m_inn = m[j] + m[i]
        pcom_inn /= m_inn
        vcom_inn /= m_inn

        orb_inn = KU.cart_to_kepl(p[j] - p[i], v[j] - v[i], m[j], m[i])
        a_inn, e_inn = orb_inn[0], orb_inn[1]

        dp_sin_inn = pcom_inn - p[k]
        dv_sin_inn = vcom_inn - v[k]
        orb_out = KU.cart_to_kepl(dp_sin_inn, dv_sin_inn, m[k], m_inn)

        a_sin_inn = orb_out[0]
        vrad_sin_inn = (dp_sin_inn * dv_sin_inn).sum()
        d_sin_inn = (dp_sin_inn * dp_sin_inn).sum() ** 0.5
        sta.d_sin_inn_ratio = d_sin_inn / a_inn
        a_apo = a_inn * (1 + e_inn)
        frel = m[i] * m[j] / (a_apo * a_apo)
        ftid = 2 * m_inn * m[k] * a_apo / (d_sin_inn * d_sin_inn * d_sin_inn)
        ffact = ftid / frel

        # print(time, d_sin_inn/a_inn, a_inn, d_sin_inn, d_sin_inn/Rhill_sum, vrad_sin_inn, inva_ik, inva_jk)
        if d_sin_inn < 2 * a_inn and sta.begin_time == -1.0:
            if sta.logfile: sta.logfile.write("t={:1.3e} : [{:d},{:d}] : d_sin_inn < 2 a_inn, "
                                              "simulation begin time\n".format(time, i, j))
            sta.begin_time = time

        if ffact > 1:
            if sta.status != "interacting":
                if sta.logfile: sta.logfile.write("t={:1.3e} : [{:d},{:d}] : binary is bound, "
                                                  "but too perturbed to be an excursion\n".format(time, i, j))
                sta.status = "interacting"
        else:
            # Unperturbed, according to some criterion
            if a_sin_inn > 0:
                # Excursion, maybe
                Pinn = 2 * np.pi * (a_inn ** 3 / m_inn) ** 0.5
                sta.average_Pinn(Pinn)
                if sta.status != "excursion" and vrad_sin_inn >= 0:
                    # An excusion can only begin with the binary-single recoiling. if it's not the case, it's not
                    sta.ex_vdotr = vrad_sin_inn
                    sta.assign_binary_and_check_exchange(i, j, k)
                    sta.ex_tbegin = time
                    if sta.logfile: sta.logfile.write("t={:1.3e} : [[{:d},{:d}],{:d}] : a_sin_inn = {:e} au > 0, "
                                                      "excursion\n".format(time, *sta.ibin, sta.isin, a_sin_inn))
                    sta.status = "excursion"
                elif sta.status == "excursion" and sta.exmap:
                    # This was already classified as an excursion
                    old_vdotr = sta.ex_vdotr
                    if vrad_sin_inn < 0.0 and old_vdotr > 0.0:
                        sta.ex_last = Excursion(m[j] + m[i], m[k],
                                                orb_inn[0], orb_inn[1], orb_inn[2], orb_inn[4],
                                                orb_out[0], orb_out[1], orb_out[2], orb_out[4],
                                                time)
                    sta.ex_vdotr = vrad_sin_inn

            else:  # a_sin_inn >= 0
                # Unbound, check for convergence
                if vrad_sin_inn >= 0:
                    # Potential escape
                    if sta.status != "escape":
                        sta.assign_binary_and_check_exchange(i, j, k)
                        if sta.logfile: sta.logfile.write(
                            "t={:1.3e} : ([{:d},{:d}],{:d}) : a_sin_inn = {:e} au < 0 and vrad > 0, "
                            "escape candidate\n".format(time, *sta.ibin, sta.isin, a_sin_inn))
                        sta.end_time = time
                        sta.status = "escape"
                    else:
                        if d_sin_inn > a_escape_mult * a_inn:
                            # It's an escape
                            # print(time, d_sin_inn / (2 * Rhill_sum), d_sin_inn/a_inn, vrad_sin_inn)
                            if sta.logfile: sta.logfile.write(
                                "t={:1.3e} : ([{:d},{:d}],{:d}) : a_sin_inn > {:g} a_inn, "
                                "escape completed\n".format(time, *sta.ibin, sta.isin, a_escape_mult))
                            sta.orb_inn = orb_inn
                            sta.orb_out = orb_out
                            return True
                else:  # vrad_sin_inn < 0
                    # Hyperbolic but converging
                    if sta.status != "will_interact":
                        sta.end = False
                        sta.assign_binary_and_check_exchange(i, j, k)
                        if sta.logfile: sta.logfile.write(
                            "t={:1.3e} : ([{:d},{:d}],{:d}) : a_sin_inn = {:e} au < 0 and vrad < 0, "
                            "binary-single will interact soon\n".format(time, *sta.ibin,
                                                                        sta.isin, a_sin_inn))
                        sta.status = "will_interact"

    if was_excursion and sta.status != "excursion":
        tex = time - sta.ex_tbegin
        avePinn = sta.ex_Pinn_Nmean[0]
        sta.ex_Pinn_Nmean = (-1, -1)

        lasted_longer_than_1period = tex > avePinn
        if lasted_longer_than_1period:
            sta.Nexcurs += 1
            if sta.exmap: sta.ex_list.append(sta.ex_last)
        else:
            if sta.logfile: sta.logfile.write("            : ended excursion : tex / Pinn = {:e}, "
                                              "invalidated\n".format(tex / avePinn))

    return False


def inv_semimajor_d(dp, dv, msum):
    v2 = (dv * dv).sum()
    r = (dp * dp).sum() ** 0.5
    invneg_a = 2 / r - v2 / msum
    return invneg_a, r


def total_energy(p, v, m):
    pot = kin = 0.0
    for key, values in ind_loop.items():
        ip, jp, kp = values
        kin += 0.5 * m[kp] * (v[kp] * v[kp]).sum()
        dp = p[ip] - p[jp]
        pot += m[ip] * m[jp] / (dp * dp).sum() ** 0.5
    etot = kin - pot
    return etot


def run_options(result=None):
    if result is None:
        result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result.add_argument("-S", help="save trajectories",
                        dest="save_traj", action="store_true")
    result.add_argument("-L", help="save logfile",
                        dest="logging", action="store_true")
    result.add_argument("-I", help="individual save",
                        dest="individual", action="store_true")
    result.add_argument("-T", help="print timing",
                        dest="timing", action="store_true")

    result.add_argument("-tfin_mult", help="final time, in units of time_to_encounter",
                        dest="tfin_mult", type=float, default=1e3)
    result.add_argument("-idist", help="initial distance, in units of binary semimajor axis",
                        dest="idist_mult", type=float, default=5)

    result.add_argument("-nogw", help="do not use gw",
                        dest="nogw", action="store_true")
    result.add_argument("-nocoll", help="do not use collisions",
                        dest="nocoll", action="store_true")
    result.add_argument("-exmap", help="excursion map",
                        dest="exmap", action="store_true")
    return result


def option_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result = run_options(result)
    result.add_argument("fname", help="run name",
                        type=str, default="testrun", nargs='?')

    result.add_argument("-m1", help="mass m1, in MSun",
                        dest="m1", type=float, default=10.0)
    result.add_argument("-m2", help="mass m2, in MSun",
                        dest="m2", type=float, default=10.0)
    result.add_argument("-m3", help="mass m3, in MSun",
                        dest="m3", type=float, default=10.0)

    result.add_argument("-ainn", help="binary semimajor axis, in au",
                        dest="a_inn", type=float, default=1.0)
    result.add_argument("-einn", help="binary eccentricity",
                        dest="e_inn", type=float, default=0.0)
    result.add_argument("-omeinn", help="binary orbit argument of pericenter, in degrees",
                        dest="ome_out", type=float, default=0.0)
    result.add_argument("-Minn", help="binary mean anomaly, in degrees",
                        dest="M_inn", type=float, default=35)

    result.add_argument("-b", help="impact parameter, in au",
                        dest="b", type=float, default=1.0)
    result.add_argument("-vinf", help="velocity at infinity, in km/s",
                        dest="vinf", type=float, default=1.0)
    result.add_argument("-iout", help="outer orbit inclination, in degrees",
                        dest="i_out", type=float, default=0.0)
    result.add_argument("-omeout", help="outer orbit argument of pericenter, in degrees",
                        dest="ome_out", type=float, default=0.0)
    result.add_argument("-Omeout", help="outer orbit longitude of ascending node, in degrees",
                        dest="Ome_out", type=float, default=0.0)

    result.add_argument("-R1", help="R1, in Rsun. Defaults to Schwarzschild radius if not provided",
                        dest="R1", type=float, default=None)
    result.add_argument("-R2", help="R2, in Rsun. Defaults to Schwarzschild radius if not provided",
                        dest="R2", type=float, default=None)
    result.add_argument("-R3", help="R3, in Rsun. Defaults to Schwarzschild radius if not provided",
                        dest="R3", type=float, default=None)
    result.add_argument("-Rschw", help="collision radius, in units of Schwarzschild radii",
                        dest="Rschw_mult", type=float, default=10)
    return result


if __name__ == "__main__":
    args = option_parser().parse_args()

    M_inn = np.radians(args.M_inn)
    nu_inn = KU.M_to_nu(M_inn, args.e_inn)

    vnb_to_kms = 2 * np.pi * KU.au2km / KU.yr2sec
    vinf = args.vinf / vnb_to_kms

    i_out = np.radians(args.i_out)
    Ome_out = np.radians(args.Ome_out) if (1e-8 < i_out < np.pi - 1e-8) else 0.0
    ome_out = np.radians(args.ome_out)
    ome_inn = np.radians(args.ome_inn) if (1e-8 < args.e_inn) else 0.0

    # Override radii, if needed
    R1 = args.Rschw_mult * 2 * args.m1 / KU.speed_of_light ** 2 if args.R1 is None else args.R1 * KU.RSun2au
    R2 = args.Rschw_mult * 2 * args.m2 / KU.speed_of_light ** 2 if args.R2 is None else args.R2 * KU.RSun2au
    R3 = args.Rschw_mult * 2 * args.m3 / KU.speed_of_light ** 2 if args.R2 is None else args.R3 * KU.RSun2au

    p, v, m, R, orb_inn, orb_out, time_to_enc = setup_system(args.m1, args.m2, args.m3,
                                                             R1, R2, R3,
                                                             args.a_inn, args.e_inn, ome_inn, nu_inn,
                                                             vinf, args.b, i_out, ome_out, Ome_out,
                                                             args.idist_mult)

    t_stop = time_to_enc * args.tfin_mult
    status = run_system(p=p, v=v, m=m, R=R,
                        ftime=t_stop, dtout=0.00001,
                        time_to_enc=time_to_enc,
                        save_traj=args.save_traj,
                        logging=args.logging,
                        fname=args.fname,
                        timing=args.timing,
                        setname=None,
                        collisions=not args.nocoll,
                        gw=not args.nogw,
                        exmap=args.exmap)

    print("\nFinal state:")
    print(status)
    if status.exmap:
        if len(status.ex_list) != status.Nexcurs:
            print("WARNING: found inconsistent number of excursion (Nex={:d}!={:d}=len(exlist)".format(status.Nexcurs,
                                                                                                       len(status.ex_list)))
        for ix, ex in enumerate(status.ex_list):
            if args.timing:
                print("Excursion #{:d}".format(ix))
                print(ex)
        exmapdb = pd.DataFrame([ex.as_dict() for ex in status.ex_list])
        exmapdb.to_hdf(args.fname + "_exmap.h5", "exmap")
