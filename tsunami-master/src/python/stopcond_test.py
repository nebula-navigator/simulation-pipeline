import tsunami
import numpy as np
import matplotlib.pyplot as plt



# Tsunami() defaults to L=au, M=MSun
code = tsunami.Tsunami()

print("Time units in yr:", code.Tscale)
print("Length units in pc:", code.Lscale)
print("Mass units in MSun:", code.Mscale)
print("Velocity units in pc/yr:", code.Vscale)

# Initial conditions - pythagorean problem
m = np.array([1., 1., 0.1])
print("Masses", m)

p = np.array([[1.,0.,0.], [-1.,0.,0.], [7.,0.,0.]])
print("Positions", p)

v = np.array([[0.,0.5,0.], [0.,-0.5,0.], [0.,0.5001,0.]])
print("Velocities", v)

r = np.array([0., 0., 0.])
print("Radii", r)

# Particle type: ignore, just leave it to -1 for every particle
st = np.array([-1, -1, -1])

# Add particle data to code. Coordinates will be rescaled to center of mass.
code.add_particle_set(p, v, m, r, st)

# Synchronize internal code coordinates with interface
code.sync_position_and_velocities(p, v)
totp = [p.copy()]

# Initialize stopping condition
code.initialize_stopping_condition(1)

# Output timestep
dt = 2

# Final time
ft = 1000000

# Current time
realt = 0
while realt < ft:
    realt = realt + dt

    # Evolve system to realt - NOTE: the system final time won't be exacly realt, but close 
    code.evolve_system(realt)

    # Synchronize realt to the real internal system time
    realt = code.ctime

    # Check if stopcond happened
    if code.stopcond:
        print("\nStopped at time", realt)
        print("Breakup occurred at", code.StopLogger.breakup_time)
        break

    # Synchronize coordinates to interface
    code.sync_position_and_velocities(p, v)
    
    totp.append(p.copy())
    print("time={:g}/{:g} - {:2.1f}% - dE/E0={:e} - E={:g} ".format(realt, ft, realt/ft*100,
                                                                    code.deltaE, code.energy),
                                                                    end="\r")

print("Instability time at ", code.StopLogger.instability_time)
print("Escaped particle index: ", code.StopLogger.escape_id)
print("Final binary")
print("         indices       :", code.StopLogger.FinalBinary.i, code.StopLogger.FinalBinary.j)
print("         semimajor axis:", code.StopLogger.FinalBinary.a)
print("         eccentricity  :", code.StopLogger.FinalBinary.e)

totp = np.vstack(totp)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.plot(totp[::3,0], totp[::3,1], lw=2.5)
ax.plot(totp[1::3,0], totp[1::3,1], lw=2.5)
ax.plot(totp[2::3,0], totp[2::3,1], lw=2.5)
plt.show()
