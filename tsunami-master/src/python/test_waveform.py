import numpy as np
import tsunami
import matplotlib.pyplot as plt

KU = tsunami.KeplerUtils()


def setup_binary(a=0.01, e=0.0, m1=20.0, m2=30.0, nu=np.pi,
                 i=0.0, ome=0.0, Ome=0.0,
                 pn1=True, pn2=True, wPNs=True):
    pos_vel2 = np.array([0., 0., 0., 0., 0., 0.])

    pos_vel1 = KU.kepl_to_cart(pos_vel2, m1, m2, a, e, i, ome, Ome, nu, pn1=pn1 and wPNs, pn2=pn2 and wPNs)

    m = np.array([m1, m2])
    p = np.array([pos_vel1[:3], pos_vel2[:3]])
    v = np.array([pos_vel1[3:], pos_vel2[3:]])

    return m, p, v


def run_binary(m, p, v, ft_P, dt_P,
               D=8e3 * KU.pc2au,  # 8 kpc, distance from source
               theta=0.0, phi=0.0,  # observing angle from source in spherical coordinates
               pn1=True, pn2=True, pn25=True, pn3=True, pn35=True, wPNs=True):

    code = tsunami.Tsunami()  # units: G=1, L=au, M=MSun
    code.Conf.wPNs = wPNs
    code.Conf.pn1 = pn1
    code.Conf.pn2 = pn2
    code.Conf.pn25 = pn25
    code.Conf.pn3 = pn3
    code.Conf.pn35 = pn35

    # 10 Schwarzschild radii as collision radii
    R = 10 * 2 * m / code.speed_of_light ** 2
    st = np.array([0, 0])

    code.add_particle_set(p, v, m, R, st)

    t = 0.0
    dp = p[0] - p[1]
    dv = v[0] - v[1]
    orb = KU.cart_to_kepl(dp, dv, m[0], m[1], pn1=pn1 and wPNs, pn2=pn2 and wPNs)
    P = 2 * np.pi * (orb[0] ** 3 / (m[0] + m[1])) ** 0.5

    hplus, hcross = code.gravitational_waveform(0, 1, theta, phi, D)

    # time, a, e, i, ome, Ome, hplus, hcross
    data = [(t * code.Tscale, orb[0], orb[1], orb[2], orb[3], orb[4], hplus, hcross)]

    dt = dt_P * P
    tfin = ft_P * P
    while t < tfin:
        t += dt
        code.evolve_system_dtmax(t)
        t = code.time  # synchronize to internal time

        code.sync_internal_state(p, v)

        if code.check_collision():
            print("Collision at time {:g} yr".format(t * code.Tscale))
            break

        dp = p[0] - p[1]
        dv = v[0] - v[1]
        orb = KU.cart_to_kepl(dp, dv, m[0], m[1], pn1=pn1 and wPNs, pn2=pn2 and wPNs)
        P = 2 * np.pi * (orb[0] ** 3 / (m[0] + m[1])) ** 0.5
        dt = dt_P * P
        hplus, hcross = code.gravitational_waveform(0, 1, theta, phi, D)
        data.append((t * code.Tscale, orb[0], orb[1], orb[2], orb[3], orb[4], hplus, hcross))

        print("Time: {:4.4f} yr, t/tfin {:2.1%}. Delta Err = {:1.4e}".format(code.time * code.Tscale, code.time / tfin,
                                                                             code.deltaE), end="\033[K\r")

    return np.array(data)


def plot(data):
    ################ Plotting ################
    f, ax = plt.subplots(3, 1, figsize=(5, 5), sharex=True, gridspec_kw=dict(hspace=0), tight_layout=True)

    t = data[:, 0] * KU.yr2sec / 60
    a = data[:, 1]
    e = data[:, 2]
    i = data[:, 3]
    ome = data[:, 4]
    Ome = data[:, 5]
    ome_un = np.unwrap(ome)

    hplus = data[:, 6]
    hcross = data[:, 7]

    ax[0].plot(t, a)
    ax[0].set_ylabel('$a$ [au]')

    ax[1].plot(t, e)
    ax[1].set_ylabel('$e$')

    ax[2].plot(t, ome_un)
    ax[2].set_ylabel(r'$\omega$')
    ax[2].set_xlabel('time [min]')

    f, ax = plt.subplots(1, 1, figsize=(5, 4), sharex=True, gridspec_kw=dict(hspace=0), tight_layout=True)
    ax.plot(t, hplus, label=r'$h^{+}$')
    ax.plot(t, hcross, label=r'$h^{\times}$')
    ax.set_ylabel("strain")
    ax.set_xlabel('time [min]')
    ax.legend()

    plt.show()

if __name__ == '__main__':
    wPNs = True
    m, p, v = setup_binary(a=0.01, e=0.0, m1=20.0, m2=30.0, nu=np.pi, i=0.0, ome=0.0, Ome=0.0, wPNs=wPNs)
    data = run_binary(m, p, v, ft_P=5, dt_P=1e-3, D=8e3 * KU.pc2au, theta=0.0, phi=0.0, wPNs=wPNs)
    plot(data)
