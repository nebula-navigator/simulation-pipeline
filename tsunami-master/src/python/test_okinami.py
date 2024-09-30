import tsunami
from tsunami import stateid as okid
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# The code uses N-body units, where G=1 and the other units are set by the equation
# Time^2 = Length^3 / (G * Mass)
# The constructor below accepts the Mass unit in Solar masses, and the Length unit in parsecs
secularcode = tsunami.Okinami(1, 1)  # Use 1 Msun, 1 au units --> 1 yr = 2pi code units
# The above is default, so secularcode = Okinami() produces same results

print("Mass unit in solar masses: {:g} MSun".format(secularcode.Mscale))
print("Length unit in parsecs: {:g} pc".format(secularcode.Lscale))
print("Time unit in years: {:g} yr\n".format(secularcode.Tscale))

np.set_printoptions(precision=16)
mpl.rcParams['figure.raise_window'] = False

# Using KeplerUtils to convert between units, convert Keplerian to Cartesian coordinates and viceversa
KU = tsunami.KeplerUtils(1, 1)

RSun = KU.RSun2au
MJup = KU.MJup2MSun
RJup = KU.RJup2au
sec2yr = 1/KU.yr2sec

# Set inner and outer semimajor axis
a1, a2 = 6, 100

# Set masses
m1, m2, m3 = 1 * MJup, 1, 40 *MJup

# Set inner and outer eccentricity
e1, e2 = 0.001, 0.6

# Set argument of pericenter
ome1, ome2 = np.radians(45), np.radians(0)

# Set inclination with respect to total angular momentum
i1, i2 = np.radians(64.7), np.radians(0.3)

# Set longitude of the ascending node
# NOTE: the equations are derived in Delaunay coordinates in the invariable plane,
#       which means that Ome1 - Ome2 = pi
# We only need to set and evolve one
Ome1, Ome2 = 0.0, np.pi

# Change to Delaunay notation for familiarity
h1 = Ome1
g1, g2 = ome1, ome2

# Compute mutual inclination and total angular momentum H from (possibly) non Delaunay coordinates
# Here Ome1 and Ome2 don't have to satisfy Ome1 - Ome2 = pi
i_mut, inc1, inc2, H = secularcode.imutual_delaunay(m1, m2, m3,
                                                    a1, a2, e1, e2,
                                                    i1, i2, Ome1, Ome2)

print("Total inclination", np.degrees(i_mut),
      "Inclination1", np.degrees(inc1),
      "inclination2", np.degrees(inc2), "\n")

# Initialize integration constants (not really constant with external derivatives)
secularcode.initialize_constants(m1, m2, a2, m3)

# Set error tolerance (default absolute_error = 1e-8, relative_error = 1e-7)
abs_err = 1e-11
rel_err = 1e-11
secularcode.set_integrator_tolerance(abs_err=abs_err, rel_err=rel_err)

# Check for stability using stability criteria
# 0 or False: no stability check
# 1 or True: use Mardling & Aarseth 2001 criteria
# 2: use Petrovich 2015 criteria (not implemented)
# 3: use AMD Hill stability criteria (Petit & Laskar 2018)
secularcode.check_stability = 0

# Enable check for collisions in the inner orbit, when a1 * (1 - e1) < R1 + R2
secularcode.check_collisions = False

# Set tidal parameters of inner orbit (body 1 and 2)
# For the tidal model, see Hut 1981 and Fabrycky & Tremaine 2007
k1 = 0.14 # Apsidal motion constant (dimensionless)
tau1 = 0.1 * sec2yr * secularcode.Tscale # Time lag (time units)
k2 = 0.3 # Apsidal motion constant (dimensionless)
tau2 = 100. * sec2yr * secularcode.Tscale  # Time lag (time units)
secularcode.set_tidal_parameters(k1, tau1, k2, tau2)

# Enable / Disable tides
secularcode.tides = False

# Enable / Disable gr (precession + decay, inner orbit only)
secularcode.gr = False

# Enable / Disable octupole term
secularcode.octupole = True

# Set radii - needed for tides and to check for collisions
R1, R2, R3 = 1,  1*RJup, 0.95*RJup

# Setup initial state array
#  0,  1,  2,  3,  4,  5, 6,  7,  8,  9, 10, 11, 12, 13
# e1, e2, g1, g2, h1, a1, H, a2, m1, m2, m3, R1, R2, R3
# Names correspond to usual Delaunay variable, where 1 and 2 refers to the inner and outer orbits
# H is the total angular momentum. See Naoz et al. 2013a for more details
y = np.array([e1, e2, g1, g2, h1, a1, H, a2, m1, m2, m3, R1, R2, R3])

# For ease of access, use the stateid enum from tsunami.py
print("State array members:")
print(list(okid), "\n")

# If there are no external derivatives injected, the following quantities remain constant:
# R1, R2, R3, m1, m2, m3, a1, a2, H
# General relativity will change: e1, a1, g1, H
# Tides will change: e1, a1, g1

# Test compute_derivatives
dydt = secularcode.compute_derivatives(y)
#print(y)
#print(dydt)

print("Now running secular evolution with Okinami\n")

time = 0.0
ftime = 2.5e7 / secularcode.Tscale  # Convert to nbody units
dtout = 1e4 / secularcode.Tscale
i_mut_l, i1_l, e1_l, t_l = [], [], [], []
while time < ftime:
    time = time + dtout
    secularcode.evolve_system(y, time)
    time = secularcode.ctime

    # Calculate i_mut from state array
    i_mut_t = secularcode.itot_from_y(y)

    # Calculate i1, i2 from i_mut
    i1, i2 = tsunami.Okinami.i1_i2_from_itot(y[okid.m1], y[okid.m2], y[okid.m3],
                                     y[okid.a1], y[okid.a2],
                                     y[okid.e1], y[okid.e2],
                                     i_mut)

    i1_l.append(np.degrees(i1))
    i_mut_l.append(np.degrees(i_mut_t))
    e1_l.append(y[okid.e1])
    t_l.append(time * secularcode.Tscale)
    print("Time (secular) = {:3.3e} yr".format(time * secularcode.Tscale), end="\033[K\r")

print("\n1-e1 maximum = {:e}, e1 minimum = {:g}\n".format(1-secularcode.e1max, secularcode.e1min))
e1_l = np.array(e1_l)

hfont = { 'fontsize' : 19, 'labelpad' : 10 }

f, ax = plt.subplots(2, sharex=True, figsize=(10, 7), gridspec_kw=dict(hspace=0),  tight_layout=True)
ax[0].plot(t_l, 1 - e1_l, lw=2.5)
ax[0].axhline(y=1 - secularcode.e1max, lw=1.5, c="black", ls=":", label="$e_{1,\\rm max}$")
ax[0].set_yscale('log')
ax[0].set_ylabel("$1-e_{1}$", **hfont)
ax[0].tick_params(labelsize=19)
ax[0].legend()

ax[1].plot(t_l, i_mut_l, label='OKINAMI (secular)', lw=2.5)
ax[1].set_ylabel("$i_{\\rm mut}$ [deg]", **hfont)
ax[1].set_xlabel("time [yr]", **hfont)
ax[1].tick_params(labelsize=19)
ax[1].legend()

plt.show(block=False)

print("Now setting up the N-body simulation with Tsunami\nGenerating cartesian coordinates for the triple system\n")

# Direct N-body requires setting true anomalies of the triple
nu1, nu2 = np.pi, 0.6*np.pi

# Setting up positions+velocities vector of body 2
pos_vel2 = np.array([0.,0.,0., 0.,0.,0.])
minn = m1+m2

# Setting up positions+velocities vector of body 1, inner orbit
pos_vel1 = KU.kepl_to_cart(pos_vel2, m1, m2, a1, e1, i1, ome1, Ome1, nu1)

# Rescale inner binary to the center of mass
pos_vel_inn_com = (pos_vel1*m1 + pos_vel2*m2)/minn
pos_vel1 -= pos_vel_inn_com
pos_vel2 -= pos_vel_inn_com

# Setup outer orbit
mout = minn+m3
pos_vel3 = KU.kepl_to_cart(np.array([0., 0., 0., 0., 0., 0.]), minn, m3, a2, e2, i2, ome2, Ome2, nu2)

# Rescale everything to the center of mass
triple_com = (pos_vel1*m1 + pos_vel2*m2 + pos_vel3*m3)/mout
pos_vel1 -= triple_com
pos_vel2 -= triple_com
pos_vel3 -= triple_com

# Initialize code with N-body units L=au, M=MSun
directcode = tsunami.Tsunami(1, 1)

# Set masses
m = np.array([m1, m2, m3])
print("m", m)

# Set positions
p = np.array([pos_vel1[:3], pos_vel2[:3], pos_vel3[:3]])
print("p", p)

# Set velocities
v = np.array([pos_vel1[3:], pos_vel2[3:], pos_vel3[3:]])
print("v", v)

# Set particle radii, needed only for collisions and tides
r = np.array([0., 0., 0.])

# Particle id, can be ignored
st = np.array([-1, -1, -1])

# Adding particle set
directcode.add_particle_set(p, v, m, r, st)

time = 0.0
ftime = 2.5e7 / directcode.Tscale
dtout = 1e4 / directcode.Tscale
i1_l, e1_l, t_l = [], [], []

plot_e = ax[0].plot(t_l, e1_l, lw=2.5)[0]
plot_i = ax[1].plot(t_l, i1_l, lw=2.5, label='TSUNAMI (direct)')[0]

input("\nPress enter to start the 3-body simulation\nIt will take about 5 minutes\n")

while time < ftime:
    time = time + dtout

    directcode.evolve_system(time)
    time = directcode.time
    directcode.sync_internal_state(p, v)
    t_l.append(time * directcode.Tscale)

    p_inn = p[0]-p[1]
    v_inn = v[0]-v[1]
    orb1 = KU.cart_to_kepl(p_inn, v_inn, m1, m2)
    p_com_inn = (p[0]*m1 + p[1]*m2)/minn
    v_com_inn = (v[0]*m1 + v[1]*m2)/minn
    p_out = p[2]-p_com_inn
    v_out = v[2]-v_com_inn

    linn = np.cross(p_inn, v_inn)
    linn = linn/(linn*linn).sum()**0.5
    lout = np.cross(p_out, v_out)
    lout = lout/(lout*lout).sum()**0.5

    i_mut = np.arccos((lout*linn).sum())
    i_mut = np.degrees(i_mut)

    e1_l.append(orb1[1])
    i1_l.append(i_mut)
    
    plot_e.set_data(t_l, 1-np.array(e1_l))
    plot_i.set_data(t_l, i1_l)
    plt.draw()
    plt.pause(1e-8)
        
    orb2 = KU.cart_to_kepl(p_out, v_out, minn, m3)
    
    print("Time (direct) = {:3.3e} yr".format(time * directcode.Tscale), end="\033[K\r")

plt.show()
