import tsunami
import numpy as np
import matplotlib.pyplot as plt


# Initialize code units, mandatory to set correct timescale and non-Newtonian forces
# Uses N-body (Henon) units, so implicitly G=1 and we need only Mass and Length units.
# Initialize: Tsunami(Mass unit in MSun, Length unit in AU)
# Tsunami() defaults to Tsunami(1, 1),
code = tsunami.Tsunami(1.0, 1.0)

print("Time units in yr:", code.Tscale)
print("Length units in au:", code.Lscale)
print("Mass units in MSun:", code.Mscale)
print("Velocity units in km/s:", code.Vscale)


################################ start Tsunami options

# Disable extra forces (already disabled by default)
code.Conf.wPNs = False  # PNs
code.Conf.wEqTides = False  # Equilibrium tides
code.Conf.wDynTides = False  # Dynamical tides

##################### To be implemented
#code.Conf.wExt = False  # External force
#code.Conf.wExt_vdep = False  # External force that is velocity dependent
#code.Conf.wMassEvol = False  # Mass evolution (loss or gain)

# Disable collisions
# dcoll multiples the radius of particles when checking for collisions
code.Conf.dcoll = 0.0

# If you change some options *AFTER* calling add_particle_set, you need to call
code.commit_parameters()  # Here is not needed

################################ end Tsunami options


# Initial conditions - Pythagorean problem
m = np.array([3., 4., 5.])
print("Masses", m)  

p = np.array([[1.,3.,0.], [-2.,-1.,0.], [1.,-1.,0.]])
print("Positions", p) 

v = np.array([[0.,0.,0.], [0.,0.,0.], [0.,0.,0.]])
print("Velocities", v)

r = np.array([0., 0., 0.])
print("Radii", r)

# Particle type: ignore, just leave it to -1 for every particle
st = np.array([-1, -1, -1])

# Add particle data to code. Coordinates will be rescaled to the center of mass (CoM).
code.add_particle_set(p, v, m, r, st)

# Synchronize internal code coordinates with the Python interface
code.sync_internal_state(p, v)  # Not really needed, but useful if provided coordinates are not CoM centered
totp = [p.copy()]

# Output timestep
dt = 0.1

# Final time
ft = 65

# Current time
time = 0
while time < ft:
    time = time + dt

    # Evolve system to realt - NOTE: the system final time won't be exacly realt, but close 
    code.evolve_system(time)

    # Synchronize realt to the real internal system time
    time = code.time

    # Returns the acceleration vectors of a pair of particles onto each other
    #a1, a2 = code.get_accelerations_of_particle_pair(0, 1)

    # Synchronize coordinates to Python interface
    code.sync_internal_state(p, v)
    
    totp.append(p.copy())
    print("time={:g}/{:g} - {:2.1%} - dE/E0={:e} - E={:g}".format(time, ft, time/ft,
                                                                  code.deltaE, code.energy),
          end="\033[K\r")

totp = np.vstack(totp)

if code.Conf.useProfiling:
    code.print_profiling()

fig = plt.figure(figsize=(5,8), tight_layout=True)
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.plot(totp[::3,0], totp[::3,1], lw=2.5)
ax.plot(totp[1::3,0], totp[1::3,1], lw=2.5)
ax.plot(totp[2::3,0], totp[2::3,1], lw=2.5)
plt.show()
