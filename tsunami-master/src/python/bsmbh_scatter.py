import matplotlib.pyplot as plt
import tsunami
import numpy as np
import argparse

keplutils = tsunami.KeplerUtils()
# Using L=au, M=MSun, G=1 (T = 1/2pi yr)

def setup_system(m1=1e6, m2=1e5, a1=100, e1=0.0, ome1=0.0, M1=0.0, Rschw_mult=1,
                 m3=10, rp2=30, i2=0.0, ome2=0.0, Ome2=0.0, d2_over_a1=50,
                 R3=0.046):
    """
    :param m1: mass SMBH 1
    :param m2: mass SMBH 2
    :param a1: semimajor axis SMBH binary
    :param e1: eccentricity SMBH binary
    :param ome1: argument of pericenter SMBH binary
    :param M1: mean anomaly SMBH binary
    :param Rschw_mult: collision radius for the SMBHs, in units of Schwartzchild radii
    :param m3: mass of the star
    :param rp2: pericenter distance of the star's parabolic orbit
    :param i2: inclination of the star's parabolic orbit
    :param ome2: argument of pericenter of the star's  parabolic orbit
    :param Ome2: longitude of the ascending node of the star's parabolic orbit
    :param d2_over_a1:
    :param R3:
    :return:
    """

    # Setup BSMBH, on x-y plane
    i1 = Ome1 = 0.0
    pos_vel2 = np.array([0.,0.,0., 0.,0.,0.])
    nu1 = keplutils.M_to_nu(M1, e1)
    minn = m1+m2
    pos_vel1 = keplutils.kepl_to_cart(pos_vel2, m1, m2, a1, e1, i1, ome1, Ome1, nu1)
    pos_vel_inn_com = (pos_vel1*m1 + pos_vel2*m2)/minn

    # Setup star on parabolic orbit
    one_minus_e2 = -1e-9  # Good enough for pericenter accuracy to 1e-8
    a2 = rp2/one_minus_e2
    e2 = 1 - one_minus_e2

    d2 = d2_over_a1 * a1
    latum = 2 * rp2
    #latum2 = a2*(1-e2*e2)
    #print(latum)
    #print(latum2)

    nu2 = latum / d2 - 1
    nu2 = nu2 / e2
    nu2 = - np.arccos(nu2)

    mout = m1+m2+m3
    pos_vel3 = keplutils.kepl_to_cart(pos_vel_inn_com, minn, m3, a2, e2, i2, ome2, Ome2, nu2)

    pos_vel_com = (pos_vel1*m1 + pos_vel2*m2 + pos_vel3*m3)/mout

    pos_vel_list = np.array([pos_vel1, pos_vel2, pos_vel3])
    pos_vel_list -= pos_vel_com

    # Set up collision radii, Schwarchild radius for BHs, tidal radius for star (G=1)
    R1 = Rschw_mult* 2*m1 / keplutils.speed_of_light**2
    R2 = Rschw_mult* 2*m2 / keplutils.speed_of_light**2
    R3 = R3

    # Use smallest BH mass, to be conservative
    #R3tide = (m2/m3)**(1/3) * R3
    # Add here GR corrections

    m = np.array([m1, m2, m3])
    R = np.array([R1, R2, R3])
    p = pos_vel_list[:,0:3].copy(order='C')  # Ensures that vectors are contiguous
    v = pos_vel_list[:,3:6].copy(order='C')

    return p, v, m, R


def run_system(p, v, m, R, ft=1000, dt_out=0.5, stop_dist=10000, save_trajectories=True, fname="0000000"):
    code = tsunami.Tsunami()

    # If we want post-Newtonians, uncomment
    #code.wPNs = True

    # If we want to disable collisions, uncomment
    #code.dcoll = 0.0

    code.add_particle_set(p, v, m, R, np.zeros_like(m, dtype=np.int64))
    code.sync_position_and_velocities(p, v)

    code.initialize_stopping_condition(2)
    
    # Initialize index of star
    code.StopLogger.set_star_id(2)
    
    # Calculate tidal radii
    Rt = (m[2]/m[:2])**(1/3) * R[2]
    Psi =(1.47 + np.exp(( m[2]-0.669 )/0.137)) / (1. + 2.34 * np.exp( (m[2]-0.669)/0.137) )
    Rt_full = Rt * Psi
    Rt_part = 2*Rt
    
    for i in [0,1]:
        code.StopLogger.set_tidal_radii(2, i, Rt_full[i], Rt_part[i])

    totp = [p.copy()]

    time = 0.0
    ftime = ft / code.Tscale
    dtout = dt_out / code.Tscale
    while time < ftime:
        time = time + dtout

        code.evolve_system(time)
        time = code.ctime
        code.sync_position_and_velocities(p, v)

        if save_trajectories:
            totp.append(p.copy())

        if code.stopcond:
            print("\nDC or FTDE time: {:g} yr".format(time*code.Tscale))
            break

        # Check for collisions, stop evolution
        # May not really trigger with our stpping conditions
        if code.CollLogger.collflag > 0:
            print("\nCollision at time: {:g} yr".format(time*code.Tscale))
            # index of colliding particles, always in ascending order
            id1 = code.CollLogger.collind[0]
            id2 = code.CollLogger.collind[1]
            break

        if escape_condition(m, p, v, stop_dist=stop_dist):
            print("\nStar ejected at time: {:g} yr".format(time*code.Tscale))
            break

        print("Time: {:3.4f} yr".format(time*code.Tscale), end="\r")

    if save_trajectories:
        totp = np.vstack(totp)
        np.savetxt(fname+".dat", totp)

    print("Direct plunge:", code.StopLogger.is_dc)
    print("Full TDE:", code.StopLogger.is_ftde)
    print("Had partial TDE:", code.StopLogger.had_ptde)
    
# This stopping condition is checked every dt_out.
# If it's not ideal, I can implement it directly in python
def escape_condition(m, p, v, stop_dist=10000, it=2):
    # Check if particle 3 (star) has escaped
    is_bound = False
    pb_com = np.zeros_like(p[0])
    vb_com = np.zeros_like(v[0])
    for i in [0,1]:
        dp = p[i] - p[it]
        dv = v[i] - v[it]

        # Keplerian orbital elements: a, e, i, ome, Ome, nu
        orb = keplutils.cart_to_kepl(dp, dv, m[i], m[it])
        is_bound = is_bound or (orb[0] > 0)

        pb_com += p[i]*m[i]
        vb_com += v[i]*m[i]

    mb = m[0] + m[1]
    pb_com/=mb
    vb_com/=mb
    dp = pb_com-p[2]
    dv = vb_com-v[2]

    orb = keplutils.cart_to_kepl(dp, dv, mb, m[2])
    is_bound = is_bound or (orb[0] > 0)
    dist = (dp*dp).sum()**0.5

    return (not is_bound) and (dist > stop_dist)


def visualize_trajectories(fname="000000"):
    pos = np.loadtxt(fname+".dat")

    f, ax = plt.subplots(figsize=(8,8))

    ax.set_aspect("equal")

    ax.plot(pos[0:-1:3,0], pos[0:-1:3,1], label="SMBH1")
    ax.plot(pos[1:-2:3,0], pos[1:-2:3,1], label="SMBH2")
    ax.plot(pos[2:-3:3,0], pos[2:-3:3,1], label="Star")
    ax.set_xlabel("x [au]")
    ax.set_ylabel("y [au]")

    plt.show()

def option_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    result.add_argument("fname", help="file name",
                        default="0000000", nargs="?")

    result.add_argument("-m1", help="SMBH1 mass, MSun",
                        dest="m1", type=float, default=1e6)
    result.add_argument("-m2", help="SMBH2 mass, MSun",
                        dest="m2", type=float, default=1e5)
    result.add_argument("-a1", help="SMBH binary semimajor axis, au",
                          dest="a1", type=float, default=100)
    result.add_argument("-e1", help="SMBH binary eccentricity",
                        dest="e1", type=float, default=0.0)
    result.add_argument("-ome1", help="SMBH binary argument of pericenter, degrees",
                        dest="ome1", type=float, default=0.0)
    result.add_argument("-M1", help="SMBH binary mean anomaly, degrees",
                        dest="M1", type=float, default=0.0)

    result.add_argument("-m3", help="star mass, MSun",
                        dest="m3", type=float, default=10)
    result.add_argument("-R3", help="star radius, RSun",
                        dest="R3", type=float, default=10)
    result.add_argument("-rp2", help="pericenter distance of stellar orbit, au",
                        dest="rp2", type=float, default=30)
    result.add_argument("-i2", help="inclination of stellar orbit, degrees",
                        dest="i2", type=float, default=0.0)
    result.add_argument("-ome2", help="argument of pericenter of stellar orbit, degrees",
                        dest="ome2", type=float, default=0.0)
    result.add_argument("-Ome2", help="longitude of ascending node of stellar orbit, degrees",
                        dest="Ome2", type=float, default=0.0)

    result.add_argument("-ft", help="final integration time, yr",
                        dest="ftime", type=float, default=1e3)
    result.add_argument("-dt", help="output timestep, yr",
                        dest="dt_out", type=float, default=0.1)
    result.add_argument("-dstop_a1", help="stopping distance, in units of binary semimajor axis",
                        dest="stop_dist_over_a1", type=float, default=100)

    return result


if __name__ == "__main__":
    args = option_parser().parse_args()

    ome1 = np.radians(args.ome1)
    M1 = np.radians(args.M1)
    ome2 = np.radians(args.ome2)
    i2 = np.radians(args.i2)
    Ome2 = np.radians(args.Ome2)
    R3 = args.R3 * keplutils.RSun2au

    stop_dist = args.stop_dist_over_a1*args.a1

    p, v, m, R = setup_system(m1=args.m1, m2=args.m2, a1=args.a1, e1=args.e1,
                              ome1=ome1, M1=M1, Rschw_mult=1,
                              m3=args.m3, rp2=args.rp2, i2=i2, ome2=ome2, Ome2=Ome2, d2_over_a1=50,
                              R3=R3)
    run_system(p, v, m, R, ft=args.ftime, dt_out=args.dt_out, stop_dist=stop_dist, save_trajectories=True, fname=args.fname)
    visualize_trajectories(fname=args.fname)