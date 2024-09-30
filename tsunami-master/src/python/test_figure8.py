import tsunami
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import matplotlib as mpl

mpl.rcParams['figure.raise_window'] = False


# Initialize units
code = tsunami.Tsunami(1, 1)

# Initial conditions - pythagorean problem
m = np.array([1., 1., 1.])
print("m", m)

r2 = np.array([-0.97000436, 0.24308753, 0.0])
r1 = -r2
r3 = np.array([0.0, 0.0, 0.0])
v3 = np.array([0.93240737, 0.86473146, 0.0])
v1 = -v3*0.5
v2 = v1

p = np.array([r1, r2, r3])
print("p", p)

v = np.array([v1, v2, v3])
print("v", v)

r = np.array([0., 0., 0.])
st = np.array([-1, -1, -1])
code.add_particle_set(p, v, m, r, st)

code.sync_internal_state(p, v)

totp = [p.copy()]

dt = 0.05
ft = 40
time = 0
while(time < ft):
    time = time + dt
    # See run_forced_timesteps.py on how the difference between evolve_system_dtmax and evolve_system
    # TLDR: evolve_system_dtmax uses the timestep given, even if smaller than the optimized one
    # Good for visualization (no jumps in trajectories), bad for accuracy
    code.evolve_system_dtmax(time)
    time = code.time
    code.sync_internal_state(p, v)
    totp.append(p.copy())
    print("time={:g}/{:g} - {:2.1%} - dE/E0={:e} - E={:g}".format(time, ft, time/ft,
                                                                  code.deltaE, code.energy),
                                                                  end="\033[K\r")

totp = np.vstack(totp)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)
ax.set_aspect('equal')
#fig.tight_layout()
fig.subplots_adjust(top=0.95, bottom=0.05, left=0.05, right=0.95)
#ax.autoscale(enable=True, axis='both', tight=True)
#ax.autoscale_view(tight=True)

baseline = 10
timeind = range(0, totp.shape[0]//3)
print("\n")

video = False
if video:
    moviewriter = FFMpegWriter(fps=40)
    moviewriter.setup(fig=fig, outfile="fig8.mp4", dpi=150)

for headt in timeind:

    ax.cla()
    tail = max(0, headt-baseline)
    
    for ip in range(3):
        ax.scatter(totp[ip+3*headt,0], totp[ip+3*headt,1], s=150)
        ax.plot(totp[ip+3*tail:ip+3*headt:3,0], totp[ip+3*tail:ip+3*headt:3,1], lw=5.5)

    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-0.5, 0.5)
    
    if video:
        print("collecting frame {:d}/{:d}".format(headt, totp.shape[0]//3), end="\r")
        moviewriter.grab_frame()
    else:
        plt.draw()
        plt.pause(0.01)
    
if video:
     moviewriter.finish()

