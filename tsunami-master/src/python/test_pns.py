import numpy as np
import tsunami
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import os

KU = tsunami.KeplerUtils()


def adot(m1, m2, a, e, pn25):
    MA = (m1 * m2 * (m1 + m2)) / ((KU.speed_of_light ** 5.) * (a ** 3.) * ((1. - e ** 2.) ** (7. / 2.)))
    eA = 1. + (73. / 24.) * (e ** 2.) + (37. / 96.) * (e ** 4.)
    dadt25 = 0.0
    if pn25: dadt25 = -(64. / 5.) * MA * eA
    return dadt25


def edot(m1, m2, a, e, pn25):
    ME = (m1 * m2 * (m1 + m2)) / ((KU.speed_of_light ** 5.) * (a ** 4.) * ((1. - e ** 2.) ** (5. / 2.)))
    eE = 1. + (121. / 304.) * (e ** 2.)
    dedt25 = 0.0
    if pn25: dedt25 = -(304. / 15.) * e * ME * eE
    return dedt25


def omedot(m1, m2, a, e, pn1, pn2):
    eta = (m1 * m2) / ((m1 + m2) ** 2.)
    w1 = w2 = 0.0

    if pn1: w1 = (3. / (KU.speed_of_light ** 2.)) * (((m1 + m2) ** (3. / 2.)) / ((a ** (5. / 2.)) * (1. - e ** 2.)))
    if pn2: w2 = - ((3. * ((m1 + m2) ** (5. / 2.))) / (4. * (KU.speed_of_light ** 4.))) * (
                (10. + 4. * eta - (1. + 10. * eta) * (e ** 2.)) / ((a ** (7. / 2.)) * ((1. - e ** 2.) ** 2.)))
    return w1 + w2


def deriv_func(t, y, m1, m2, pn1, pn2, pn25, pn3, pn35):
    a, e, ome = y
    return [adot(m1, m2, a, e, pn25), edot(m1, m2, a, e, pn25), omedot(m1, m2, a, e, pn1, pn2)]


def setup_binary(a=0.1, e=0.99, m1=50.0, m2=50.0, nu=np.pi, pn1=True, pn2=True):
    i = ome = Ome = 0.0 # rad
    pos_vel2 = np.array([0.,0.,0., 0.,0.,0.])

    pos_vel1 = KU.kepl_to_cart(pos_vel2, m1, m2, a, e, i, ome, Ome, nu, pn1=pn1, pn2=pn2)

    m = np.array([m1, m2])
    p = np.array([pos_vel1[:3], pos_vel2[:3]])
    v = np.array([pos_vel1[3:], pos_vel2[3:]])

    return m, p, v


def run_binary(m, p, v, ft, dt_P, pn1=True, pn2=True, pn25=True, pn3=False, pn35=False):
    code = tsunami.Tsunami()  # units: G=1, L=au, M=MSun
    # Use 1PN + 2PN + 2.5PN
    code.Conf.wPNs = True
    code.Conf.pn1 = pn1
    code.Conf.pn2 = pn2
    code.Conf.pn25 = pn25
    code.Conf.pn3 = pn3
    code.Conf.pn35 = pn35

    R = 10 * 2*m / code.speed_of_light**2
    st = np.array([0, 0])

    code.add_particle_set(p, v, m, R, st)

    t = 0.0
    t_l, e_l, a_l, ome_l = [], [], [], []
    dp = p[0]-p[1]
    dv = v[0]-v[1]
    orb = KU.cart_to_kepl(dp, dv, m[0], m[1], pn1=pn1, pn2=pn2)
    P = 2*np.pi * (orb[0]**3 / (m[0]+m[1]))**0.5
    t_l.append(t * code.Tscale)
    a_l.append(orb[0])
    e_l.append(orb[1])
    ome_l.append(orb[3])

    dt = dt_P * P
    tfin = ft / code.Tscale

    print("\nEvolving binary")

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
        orb = KU.cart_to_kepl(dp, dv, m[0], m[1], pn1=pn1, pn2=pn2)
        P = 2*np.pi * (orb[0]**3 / (m[0]+m[1]))**0.5
        dt = dt_P * P

        t_l.append(t * code.Tscale)
        a_l.append(orb[0])
        e_l.append(orb[1])
        ome_l.append(orb[3])

        print("Time: {:4.4f} yr, t/tfin {:2.1%}. Delta Err = {:1.4e}".format(code.time * code.Tscale, code.time/tfin, code.deltaE), end="\033[K\r")

    return np.array(t_l), np.array(a_l), np.array(e_l), np.array(ome_l)


def do_comparison(a, e, ome, m1, m2, ft, pn1, pn2, pn25, pn3, pn35):
    print("PN corrections:\n  PN1 = {}\n  PN2 = {}\n  PN2.5 = {}\n  PN3 = {}\n  PN3.5 = {}".format(pn1, pn2, pn25, pn3, pn35))


    ################ Analytic run ################
    t_span = [0., ft / KU.Tscale]  # yr [start, stop]
    dt = np.linspace(0., ft / KU.Tscale, 200)
    y0 = [a, e, ome]
    R1, R2 = 2 * m1 / KU.speed_of_light ** 2, 2 * m2 / KU.speed_of_light ** 2
    sumR = R1+R2

    def stop_afin(t, y, m1, m2, pn1, pn2, pn25, pn3, pn35):
        return y[0]*(1-y[1]) - sumR

    stop_afin.terminal = True

    sol1 = solve_ivp(deriv_func, t_span, y0, t_eval=dt, args=(m1, m2, pn1, pn2, pn25, pn3, pn35), events=stop_afin)
    y_ana = sol1.y
    t_ana = sol1.t

    ################ N-body run ################
    txtsave = "test_pns_nbody_pn1={:d}-pn2={:d}-pn25={:d}-pn3={:d}-pn35={:d}.txt".format(pn1, pn2, pn25, pn3, pn35)
    if not os.path.isfile(txtsave):
        m, p, v = setup_binary(a=a_0, e=e_0, m1=m1, m2=m2, pn1=pn1, pn2=pn2)
        t_nb, a_nb, e_nb, ome_nb = run_binary(m, p, v, ft=ft, dt_P=1.5, pn1=pn1, pn2=pn2, pn25=pn25, pn3=pn3, pn35=pn35)
        np.savetxt(txtsave, np.array((t_nb, a_nb, e_nb, ome_nb)).T)
    else:
        t_nb, a_nb, e_nb, ome_nb = np.loadtxt(txtsave).T

    ome_nb = np.unwrap(ome_nb)
    y_ana[2] = np.degrees(y_ana[2])

    return np.vstack((t_ana * KU.Tscale, y_ana)), np.array((t_nb, a_nb, e_nb, np.degrees(ome_nb)))


################ Initial conditions ################
ft = 160.  # yr
a_0 = 0.1  # au
e_0 = 0.99
ome_0 = 0.  # radian
m1, m2 = 50., 50.  # Msun

yana_pn35, ynb_pn35 = do_comparison(a_0, e_0, ome_0, m1, m2, ft, pn1=True, pn2=True, pn25=True, pn3=True, pn35=True)
yana_pn3, ynb_pn3 = do_comparison(a_0, e_0, ome_0, m1, m2, ft, pn1=True, pn2=True, pn25=True, pn3=True, pn35=False)
yana_pn25, ynb_pn25 = do_comparison(a_0, e_0, ome_0, m1, m2, ft, pn1=True, pn2=True, pn25=True, pn3=False, pn35=False)
yana_pn2, ynb_pn2 = do_comparison(a_0, e_0, ome_0, m1, m2, ft, pn1=True, pn2=True, pn25=False, pn3=False, pn35=False)
yana_pn1, ynb_pn1 = do_comparison(a_0, e_0, ome_0, m1, m2, ft, pn1=True, pn2=False, pn25=False, pn3=False, pn35=False)


################ Plotting ################
f, ax = plt.subplots(3, 1, figsize=(10,9), sharex=True, gridspec_kw=dict(hspace=0),  tight_layout=True)

ax1 = ax[0]
#ax1.plot(yana_pn35[0], yana_pn35[1], label='Secular average 3.5PN', lw=3.0, alpha=0.66, c="tab:purple", ls="--")
ax1.plot(ynb_pn35[0], ynb_pn35[1], label='TSUNAMI 3.5PN', lw=3.0, alpha=0.66, c="tab:purple")
#ax1.plot(yana_pn3[0], yana_pn3[1], label='Secular average 3PN', lw=3.0, alpha=0.66, c="tab:red", ls="--")
ax1.plot(ynb_pn3[0], ynb_pn3[1], label='TSUNAMI 3PN', lw=3.0, alpha=0.66, c="tab:red")
ax1.plot(yana_pn25[0], yana_pn25[1], label='Secular average 2.5PN', lw=3.0, alpha=0.66, c="tab:blue", ls="--")
ax1.plot(ynb_pn25[0], ynb_pn25[1], label='TSUNAMI 2.5PN', lw=3.0, alpha=0.66, c="tab:blue")
ax1.plot(yana_pn2[0], yana_pn2[1], label='Secular average 2PN', lw=3.0, alpha=0.66, c="tab:orange", ls="--")
ax1.plot(ynb_pn2[0], ynb_pn2[1], label='TSUNAMI 2PN', lw=3.0, alpha=0.66, c="tab:orange")
ax1.plot(yana_pn1[0], yana_pn1[1], label='Secular average 1PN', lw=3.0, alpha=0.66, c="tab:green", ls="--")
ax1.plot(ynb_pn1[0], ynb_pn1[1], label='TSUNAMI 1PN', lw=3.0, alpha=0.66, c="tab:green")
ax1.set_yscale("log")
ax1.legend()
ax1.set_ylabel('$a$ [au]')

ax2 = ax[1]
#ax2.plot(yana_pn35[0], yana_pn35[2], label='Secular average 3.5PN', lw=3.0, alpha=0.66, c="tab:purple", ls="--")
ax2.plot(ynb_pn35[0], 1-ynb_pn35[2], label='TSUNAMI 3.5PN', lw=3.0, alpha=0.66, c="tab:purple")
#ax2.plot(yana_pn3[0], 1-yana_pn3[2], label='Secular average 3PN', lw=3.0, alpha=0.66, c="tab:red", ls="--")
ax2.plot(ynb_pn3[0], 1-ynb_pn3[2], label='TSUNAMI 3PN', lw=3.0, alpha=0.66, c="tab:red")
ax2.plot(yana_pn25[0], 1-yana_pn25[2], label='Secular average 2.5PN', lw=3.0, alpha=0.66, c="tab:blue", ls="--")
ax2.plot(ynb_pn25[0], 1-ynb_pn25[2], label='TSUNAMI 2.5PN', lw=3.0, alpha=0.66, c="tab:blue")
ax2.plot(yana_pn2[0], 1-yana_pn2[2], label='Secular average 2PN', lw=3.0, alpha=0.66, c="tab:orange", ls="--")
ax2.plot(ynb_pn2[0], 1-ynb_pn2[2], label='TSUNAMI 2PN', lw=3.0, alpha=0.66, c="tab:orange")
ax2.plot(yana_pn1[0], 1-yana_pn1[2], label='Secular average 1PN', lw=3.0, alpha=0.66, c="tab:green", ls="--")
ax2.plot(ynb_pn1[0], 1-ynb_pn1[2], label='TSUNAMI 1PN', lw=3.0, alpha=0.66, c="tab:green")
#ax2.set_yscale('log')
ax2.set_ylabel('$1 - e$')

ax3 = ax[2]
#ax3.plot(yana_pn35[0], yana_pn35[3], label='Secular average 3.5PN', lw=3.0, alpha=0.66, c="tab:purple", ls="--")
ax3.plot(ynb_pn35[0], ynb_pn35[3], label='TSUNAMI 3.5PN', lw=3.0, alpha=0.66, c="tab:purple")
#ax3.plot(yana_pn3[0], yana_pn3[3], label='Secular average 3PN', lw=3.0, alpha=0.66, c="tab:red", ls="--")
ax3.plot(ynb_pn3[0], ynb_pn3[3], label='TSUNAMI 3PN', lw=3.0, alpha=0.66, c="tab:red")
ax3.plot(yana_pn25[0], yana_pn25[3], label='Secular average 2.5PN', lw=3.0, alpha=0.66, c="tab:blue", ls="--")
ax3.plot(ynb_pn25[0], ynb_pn25[3], label='TSUNAMI 2.5PN', lw=3.0, alpha=0.66, c="tab:blue")
ax3.plot(yana_pn2[0], yana_pn2[3], label='Secular average 2PN', lw=3.0, alpha=0.66, c="tab:orange", ls="--")
ax3.plot(ynb_pn2[0], ynb_pn2[3], label='TSUNAMI 2PN', lw=3.0, alpha=0.66, c="tab:orange")
ax3.plot(yana_pn1[0], yana_pn1[3], label='Secular average 1PN', lw=3.0, alpha=0.66, c="tab:green", ls="--")
ax3.plot(ynb_pn1[0], ynb_pn1[3], label='TSUNAMI 1PN', lw=3.0, alpha=0.66, c="tab:green")
#ax3.set_yscale("log")
ax3.set_ylabel('$\omega$ [deg]')

for axx in ax: axx.grid()

plt.show()
