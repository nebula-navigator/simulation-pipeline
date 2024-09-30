import tsunami
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.lines as mlines

# Initialize units
# code = tsunami.Tsunami(1, 1) # Mass (in MSun), Length (in au).
code = tsunami.Tsunami()  # If empty, defaults to MSUN-AU

# Enable postnewtonians
code.Conf.wPNs = True

# Initial conditions - cold 4-body problem scattering
mass = np.array([10., 15., 20., 25.])
print("mass", mass)

pos = np.array([[0.12155378, 3, 0.],
                [-2, 2, 0.],
                [2, -2, 0.],
                [-4, -4, 0.]])
print("position", pos)

# Particles at rest
vel = np.array([[0., 0., 0.],
                [0., 0., 0.],
                [0., 0., 0.],
                [0., 0., 0.]])
print("velocity", vel)

print("c = {:e} au/yr".format(code.speed_of_light / code.Tscale))
# G is 1
# 10x Schwarzschild radius
rad = 10 * 2 * mass / code.speed_of_light ** 2
# rad = np.array([0., 0., 0., 0.])

print("radius", rad)

# Particle id, let's use it to identify our particles. Must be integer
st = np.array([0, 1, 2, 3])

# Adding particles
code.add_particle_set(pos, vel, mass, rad, st)

# Synchronizing to the internal positions and velocities (rescaled to COM)
code.sync_internal_state(pos, vel)

totpos = [pos.copy()]
totpos_old = None

dt_out = 1e-3 / code.Tscale  # Diagnostic timestep
ft = 1.3 / code.Tscale  # Final time

coll = False
time = 0.0
while time < ft:
    time = time + dt_out

    code.evolve_system(time)
    time = code.time  # synchronize to internal time
    print("Time: {:3.4f} yr".format(time * code.Tscale), end="\033[K\r")

    code.sync_internal_state(pos, vel)

    totpos.append(pos.copy())

    if code.check_collision() > 0:
        print("\nCollision detected")
        # Collision was detected, the code has rollbacked to the last timestep before the collision
        coll = True
        id1, id2 = code.get_collision_indices()
        t_coll = code.get_collision_time()
        # You can use 
        # code.reset_collision_status()
        # to reset the collision status, and resume integration from the last timestep
        print("Collision time:", t_coll)
        print("Rollback time: ", code.time)

        mass_merg = mass[id1] + mass[id2]
        pos_merg = (pos[id1] * mass[id1] + pos[id2] * mass[id2]) / mass_merg
        vel_merg = (vel[id1] * mass[id1] + vel[id2] * mass[id2]) / mass_merg

        # we are going to remove a particle
        # code need to be reset
        # let's make a new instance
        del code
        code = tsunami.Tsunami()
        code.Conf.wPNs = True

        # delete old particles
        pos = np.delete(pos, [id1, id2], axis=0)
        vel = np.delete(vel, [id1, id2], axis=0)
        mass = np.delete(mass, [id1, id2])
        rad = np.delete(rad, [id1, id2])

        oldst = st.copy()
        st = np.delete(st, [id1, id2])

        # add new particle
        pos = np.append(pos, pos_merg.reshape((1, -1)), axis=0)
        vel = np.append(vel, vel_merg.reshape((1, -1)), axis=0)
        mass = np.append(mass, mass_merg)
        rad = np.append(rad, 10 * 2 * mass_merg / code.speed_of_light ** 2)
        mergerst = -1
        st = np.append(st, mergerst)

        print("mass", mass)
        print("position", pos)
        print("velocity", vel)

        code.add_particle_set(pos, vel, mass, rad, st)
        id_collproduct = code.N - 1  # Merger product is the last particle

        code.sync_internal_state(pos, vel)

        # set time to old one
        code.time = time

        totpos_old = totpos.copy()
        totpos = [pos.copy()]

if not coll:
    totpos_old = totpos
    oldst = st

totpos = np.vstack(totpos)
totpos_old = np.vstack(totpos_old)

colordict = {0: "tab:blue", 1: "tab:orange", 2: "tab:green", 3: "tab:purple"}

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111)
ax.set_aspect('equal')

for ip, stype in enumerate(oldst):
    ax.plot(totpos_old[ip::4, 0], totpos_old[ip::4, 1], lw=2.5, c=colordict[stype])
    ax.scatter(totpos_old[ip, 0], totpos_old[ip, 1], marker="D", s=100, c=colordict[stype])
    if not coll:
        ntimes = totpos_old.shape[0] - 4
        ax.scatter(totpos_old[ntimes + ip, 0], totpos_old[ntimes + ip, 1], marker="o", s=100, c=colordict[stype])

start = mlines.Line2D([], [], color='black', marker='D', linestyle='None', markersize=10, label='Begin')
end = mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, label='End')
marks = [start, end]

if coll:
    clmark = ax.scatter(totpos[id_collproduct, 0], totpos[id_collproduct, 1], marker="x",
                        label="Collision", s=200, c="black")
    marks.append(clmark)
    merger = mlines.Line2D([], [], color='tab:red', marker='none', linestyle='-', lw=2.5, label='Merger product')
    marks.append(merger)

    colordict[mergerst] = "tab:red"
    for ip, stype in enumerate(st):
        ax.plot(totpos[ip::3, 0], totpos[ip::3, 1], lw=2.5, c=colordict[stype])
        ntimes = totpos.shape[0] - 3
        ax.scatter(totpos[ntimes + ip, 0], totpos[ntimes + ip, 1], marker="o", s=100, c=colordict[stype])

ax.legend(handles=marks, loc="best")

plt.tight_layout()
plt.show()
