import tsunami
import numpy as np
import matplotlib.pyplot as plt


def setup_system(Mstar, Rstar, planets_list):
    """
    :param Mstar: mass of star
    :param Rstar: radius of star
    :param planets_list: list of masses + radius + orbital parameters (m, R, a, e, i, ome, Ome, nu) for each planet.
                           must be sorted by semimajor axis
    :return: p, v, m:  position, velocity and mass vector
    """
    KU = tsunami.KeplerUtils()
    m_list = [Mstar]
    R_list = [Rstar]

    # position-velocity array
    pv_star = np.array([0, 0, 0, 0, 0, 0])
    pv_com = pv_star.copy()

    pv_list = [pv_star]
    menclosed = Mstar
    print("Planets_list:\n#   {:8s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s} {:8s}".format("m [MSun]", "R [au]", "a [au]", "e",
                                                                                    "i [rad]", "ome [rad]", "Ome [rad]",
                                                                                    "M [rad]"))
    for ip, planet_prop in enumerate(planets_list):
        print("#{:d} {:8.3g} {:8.3g} {:8.3g} {:8.3g} {:8.3g} {:8.3g} {:8.3g} {:8.3g}".format(ip, *planet_prop))

        m, R, a, e, i, ome, Ome, M = planet_prop
        nu = KU.M_to_nu(M, e)

        mtot = menclosed + m

        # Add new particle around center of mass of the previous
        pos_vel_p = KU.kepl_to_cart(pv_com, menclosed, m, a, e, i, ome, Ome, nu)

        pv_list.append(pos_vel_p.copy())
        m_list.append(m)
        R_list.append(R)

        pv_com = (pos_vel_p * m + pv_com * menclosed) / mtot
        menclosed = mtot

    # Center to COM
    pv_list = np.array(pv_list)
    pv_list -= pv_com

    m = np.array(m_list)
    R = np.array(R_list)
    p = pv_list[:, 0:3].copy(order='C')  # Ensures that vectors are contiguous
    v = pv_list[:, 3:6].copy(order='C')

    return m, R, p, v


def calculate_orbits(m, p, v):
    orbit_params = []
    for ip in range(1, len(m)):
        dp, dv = p[ip] - p[0], v[ip] - v[0]
        orbpars = KU.cart_to_kepl(dp, dv, m[ip], m[0])  # Returns a, e, i, ome, Ome, nu
        orbit_params.append(orbpars)
    return orbit_params


def evolve_planetary_system(Mstar, Rstar, planet_list, dtout, ftime, mute=False):
    m, R, p, v = setup_system(Mstar, Rstar, planet_list)

    code = tsunami.Tsunami()  # Nbody units: MSun, au
    code.add_particle_set(p, v, m, R, np.zeros_like(m, dtype=np.int64))
    code.sync_internal_state(p, v)

    time = 0.0
    totp = [p.copy()]
    orbit_params = calculate_orbits(m, p, v)
    orbit_params.insert(0, [time])
    orbparams_l = [np.concatenate(orbit_params)]

    dt = dtout / code.Tscale
    final_time = ftime / code.Tscale
    while time < final_time:
        time = time + dt

        code.evolve_system(time)
        time = code.time

        code.sync_internal_state(p, v)

        totp.append(p.copy())
        orbit_params = calculate_orbits(m, p, v)
        orbit_params.insert(0, [time * code.Tscale])
        orbparams_l.append(np.concatenate(orbit_params))

        if code.check_collision() > 0:
            # index of colliding particles, always in ascending order
            id1, id2 = code.get_collision_indices()
            print(f"\nCollision detected between planet #{id1} and #{id2}")
            break

        if not mute: print(
            "time={:4.3e}/{:4.3e} yr - {:2.1%} - dE/E0={:e} - E={:g} ".format(time, final_time, time / final_time,
                                                                              code.deltaE, code.energy), end="\033[K\r")

    totp = np.vstack(totp)
    orbparams_l = np.vstack(orbparams_l)

    return totp, orbparams_l


def plot_evolution(orbpars):
    # get the number of planets from shape of orbpars columns:
    # time  (a, e, i, ome, Ome, nu) x Nplanets
    Np = (orbpars.shape[1] - 1) // 6

    t = orbpars[:, 0]

    fig, ax = plt.subplots(2, 1, figsize=(5, 4), sharex=True)

    # Plot semimajor axis, eccentricity
    for ip in range(Np):
        a = orbpars[:, 1 + ip * 6]
        e = orbpars[:, 2 + ip * 6]

        ax[0].plot(t, a, lw=2.5, label=f"planet {ip}")
        ax[1].plot(t, e, lw=2.5, label=f"planet {ip}")

    ax[0].set_ylabel("Semimajor axis [au]")
    ax[1].set_ylabel("Eccentricity")

    ax[1].set_xlabel("Time [yr]")
    for aaxx in ax:
        aaxx.grid(True)
        aaxx.legend()

    fig.subplots_adjust(hspace=0)

    plt.show()


def plot_trajectories(pos):
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))

    # Plot trajectories
    ax.set_aspect('equal')
    ax.scatter(pos[::3, 0], pos[::3, 1], alpha=0.5, marker=".", s=1)
    ax.scatter(pos[1::3, 0], pos[1::3, 1], alpha=0.5, marker=".", s=1)
    ax.scatter(pos[2::3, 0], pos[2::3, 1], alpha=0.5, marker=".", s=1)
    ax.set_xlabel("X [au]")
    ax.set_ylabel("Y [au]")
    ax.grid(True)

    plt.show()


if "__main__" == __name__:
    KU = tsunami.KeplerUtils()

    # Units: MSun, au, yr/2pi, G=1
    MEarth = KU.MEarth2MSun
    MJup = KU.MJup2MSun
    REarth = KU.REarth2au
    RJup = KU.RJup2au
    RSun = KU.RSun2au

    # Star mass and radius in msun au
    Mstar, Rstar = 1, RSun

    m1 = MEarth
    R1 = REarth
    a1 = 1
    e1 = 0.02
    M1 = np.pi / 2

    m2 = MJup
    R2 = RJup
    a2 = 4
    e2 = 0.8
    i2 = np.radians(70)
    ome2 = np.radians(150)
    Ome2 = np.radians(220)
    M2 = 3 * np.pi / 2

    planet1 = [m1, R1, a1, e1, 0.0, 0.0, 0.0, M1]
    planet2 = [m2, R2, a2, e2, i2, ome2, Ome2, M2]
    planet_list = [planet1, planet2]

    planet_list.sort(key=lambda x: x[2])

    # period planet 1 in yr
    P1 = (a1 ** 3 / (Mstar + m1)) ** 0.5 * KU.Tscale

    positions, orbpars = evolve_planetary_system(Mstar, Rstar, planet_list, dtout=P1 * 10, ftime=P1 * 1e5)
    plot_evolution(orbpars)
    plot_trajectories(positions)
