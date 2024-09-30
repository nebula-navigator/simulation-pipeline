"""
This scripts shows how to force a smaller timestep in TSUNAMI.
This might be required for making pretty videos with smooth trajectories
or simply have more data points

Important notes:
    * The timestep won't be exactly the same are requested, there will still be some imprecision
    * SHOULDN'T BE USED FOR PRODUCTION! Forcing a smaller timestep will increase numerical errors in the long run
"""

import tsunami
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from run_planetary_system import setup_system

KU = tsunami.KeplerUtils()
MJup = KU.MJup2MSun
RJup = KU.RJup2au
RSun = KU.RSun2au

# Setting up the system. A Sun-like star with two Jupiter-like planets. We use setup_system from run_planetary_system.py
# (m, R, a, e, i, ome, Ome, nu)
planets_list = [(MJup, RJup, 1, 0.0, 0.0, 0.0, 0.0, np.pi),
                (MJup, RJup, 5, 0.3, np.pi*0.15, np.pi/3, np.pi/2, 1.5*np.pi)]

m, R, p, v = setup_system(Mstar=1.0, Rstar=RSun, planets_list=planets_list)

# Save initial conditions, we will run it twice
p0, v0 = p.copy(), v.copy()

# Output timestep = 1 day
dt = 1 / KU.Tscale / 365.25
# Final time = 200 years
ft = 200 / KU.Tscale

# Estimated number of timesteps
Ndt_req = ft / dt
print("Number of timesteps requested:", int(Ndt_req))

# First run without forcing the timestep
code = tsunami.Tsunami()
code.add_particle_set(p0, v0, m, R, np.zeros_like(m, dtype=np.int64))
print("== Run without forced timestep ==")
code.sync_internal_state(p, v)
totp, timel = [p.copy()], [0.0]
time = 0.0
while time < ft:
    time = time + dt

    code.evolve_system(time)
    time = code.time

    # Synchronize coordinates to interface
    code.sync_internal_state(p, v)

    totp.append(p.copy())
    timel.append(time * KU.Tscale)
    print("time={:g}/{:g} - {:2.1%} - dE/E0={:e} - E={:g} - ds={:g} - dt={:g}".format(
        time, ft, time / ft, code.deltaE, code.energy, code.timestep, code.dtphysical), end="\033[K\r")

Ndt = len(totp)
print("\nNumber of timesteps done:     ", Ndt)
totp = np.vstack(totp)
timel= np.array(timel)


# Second run with forcing the timestep
code = tsunami.Tsunami()
code.add_particle_set(p0, v0, m, R, np.zeros_like(m, dtype=np.int64))
print("== Run with forced timestep ==")
code.sync_internal_state(p, v)

totp_forced, timel_forced = [p.copy()], [0.0]
time = 0.0
while time < ft:
    # We set the next step to be very small, so we are sure we won't overshoot the requested time
    # The minimum timestep is controlled by code.timestep
    time = time + dt

    code.evolve_system_dtmax(time)
    time = code.time

    # Synchronize coordinates to interface
    code.sync_internal_state(p, v)

    totp_forced.append(p.copy())
    timel_forced.append(time * KU.Tscale)
    print("time={:g}/{:g} - {:2.1%} - dE/E0={:e} - E={:g} - ds={:g} - dt={:g}".format(
        time, ft, time / ft, code.deltaE, code.energy, code.timestep, code.dtphysical), end="\033[K\r")

Ndt = len(totp_forced)
print("\nNumber of timesteps done:     ", Ndt)
totp_forced = np.vstack(totp_forced)
timel_forced = np.array(timel_forced)

f, ax = plt.subplots(figsize=(8, 8), tight_layout=True)

tplot = 2
for ip in range(code.N):
    Ntp = np.argwhere(timel > tplot)[0][0]
    slice = np.s_[ip:Ntp*code.N+ip:code.N]
    ax.plot(totp[slice,0], totp[slice,1], lw=8, alpha=0.5)
ax.set_prop_cycle(None)

for ip in range(code.N):
    Ntp = np.argwhere(timel_forced > tplot)[0][0]
    slice = np.s_[ip:Ntp*code.N+ip:code.N]
    ax.plot(totp_forced[slice,0], totp_forced[slice,1], lw=2, ls="--")

ax.set_ylabel("Y [au]")
ax.set_xlabel("X [au]")
ax.legend([Line2D([0], [0], color="black", lw=8, alpha=0.5),
           Line2D([0], [0], color="black", lw=2, ls="--")], ["No forced timestep", "Forced timestep"])
plt.show()
