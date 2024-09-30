import tsunami
import numpy as np
import matplotlib.pyplot as plt
from equilibrium_tide_averaged import EqTides, make_tsunami_ic


def get_spinrate_and_obliquity(dp, dv, spin1_vec, spin2_vec):
    h_vec = np.cross(dp, dv)
    h = (h_vec*h_vec).sum()**0.5
    h_hat = h_vec/h

    spin1 = (spin1_vec*spin1_vec).sum()**0.5
    spin2 = (spin2_vec*spin2_vec).sum()**0.5

    spin1_hat = spin1_vec/spin1
    spin2_hat = spin2_vec/spin2

    obl1 = np.degrees(np.arccos((spin1_hat * h_hat).sum()))
    obl2 = np.degrees(np.arccos((spin2_hat * h_hat).sum()))
    return spin1, obl1, spin2, obl2


KU = tsunami.KeplerUtils()

############# Initial conditions

print("# Setting up system")

m1 = 1  # MSun
m2 = 1  # MSun
R1 = 1 * KU.RSun2au
R2 = 1 * KU.RSun2au

# Apsidal motion constant
kap = 0.05
k1 = kap
k2 = kap

# time-lag, first in seconds, converted to N-body
tausec = 1e3
tau1 = tausec / KU.yr2sec / KU.Tscale
tau2 = tausec / KU.yr2sec / KU.Tscale

# Gyration radius, given inertia I = M * (rg * R)^2
rg = 0.28
rg1 = rg
rg2 = rg

# Initial orbit
a0, e0 = 0.09, 0.3

# Spin, first in period (days) then in angular frequency
Ps1 = 100 / 365.25 / KU.Tscale
spin1 = 2 * np.pi / Ps1

Ps2 = 0.0015 / KU.Tscale
spin2 = 2 * np.pi / Ps2

############# Run with averaged equations

print("# Running with averaged equation")

ET = EqTides(m1, m2, R1, R2, k1, k2, tau1, tau2, rg1, rg2, rotdist=False, usespin=True)

ET.setup_orbital_vectors(a0, e0)

obl1 = np.radians(30)
obl2 = 0.0

ET.setup_spin_vectors(spin1, spin2, obl1, obl2)
ET.setup_vectors()
ET.initialize_integrator()

# Copy the spin values here, before the evolution
m, R, p, v, spin, pt = make_tsunami_ic(m1, m2, R1, R2, a0, e0, ET.spin1_vec, ET.spin2_vec)

# Final time, in units of initial orbital period (in N-body units)
tfin = 5e5 * ET.P0
dt = 5e1 * ET.P0

# Evolve with averaged equation
secdata = ET.evolve(tfin, dt)

# Calculate data to plot
tsec = secdata[:, 0] / ET.P0
asec = secdata[:, 1]
esec = secdata[:, 2]
n = (ET.mtot / asec ** 3) ** 0.5 / ET.KU.Tscale

spin1sec = secdata[:, 3] / ET.KU.Tscale
spin2sec = secdata[:, 4] / ET.KU.Tscale
obl1sec = np.degrees(secdata[:, 5])
obl2sec = np.degrees(secdata[:, 6])
nsync = EqTides.pseudosync_spin(esec, n)
necc = EqTides.eccentricity_spin(esec, n)

############# Running with Tsunami

print("# Running with Tsunami")

# Initialize code
code = tsunami.Tsunami()
code.Conf.wEqTides = True

print(spin)
code.add_particle_set(p, v, m, R, pt, spin)

kaps = np.array([k1, k2])
taulag = np.array([tau1, tau2])
polyind = np.array([0.0, 0.0])  # only for dynamical tide, null values
gyrad = np.array([rg1, rg2])

code.initialize_tidal_parameters(kaps, taulag, polyind, gyrad)

code.sync_internal_state(p, v, spin)
dp = p[1]-p[0]
dv = v[1]-v[0]
sp1, obli1, sp2, obli2 = get_spinrate_and_obliquity(dp, dv, spin[0], spin[1])
orb = KU.cart_to_kepl(dp, dv, m[0], m[1])

time = 0.0
spin1nb = [sp1 / code.Tscale]
spin2nb = [sp2 / code.Tscale]
obl1nb = [obli1]
obl2nb = [obli2]
anb = [orb[0]]
enb = [orb[1]]
P0 = 2*np.pi * (orb[0]*orb[0]*orb[0]/m.sum())**0.5
tnb = [time / P0]

while time < tfin:
    time = time + dt
    
    code.evolve_system(time)
    time = code.time
    code.sync_internal_state(p, v, spin)

    dp = p[1]-p[0]
    dv = v[1]-v[0]
    orb = KU.cart_to_kepl(dp, dv, m[0], m[1])
    sp1, obli1, sp2, obli2 = get_spinrate_and_obliquity(dp, dv, spin[0], spin[1])

    spin1nb.append(sp1 / code.Tscale)
    spin2nb.append(sp2 / code.Tscale)
    obl1nb.append(obli1)
    obl2nb.append(obli2)
    anb.append(orb[0])
    enb.append(orb[1])
    tnb.append(time / P0)
    print("time={:g}/{:g} - {:2.1%} - dE/E0={:e} - E={:g} ".format(time, tfin, time/tfin,
                                                                    code.deltaE, code.energy),
                                                                    end="\033[K\r")

f, ax = plt.subplots(4, 1, figsize=(7,8), sharex=True, gridspec_kw=dict(hspace=0),  tight_layout=True)
plt.rc('lines', linewidth=2.5)
plt.rc('text', usetex=True)

ax[0].plot(tnb, anb, alpha=0.66, label="$\\textsc{tsunami}$")
ax[0].plot(tsec, asec, alpha=0.66, ls="--", c="tab:blue", label="Secular")
ax[0].set_ylabel("semimajor axis [au]")
ax[0].legend()

ax[1].plot(tnb, enb, alpha=0.66)
ax[1].plot(tsec, esec, alpha=0.66, ls="--", c="tab:blue")
ax[1].set_ylabel("eccentricity")
ax[1].set_xlabel("t [yr/2pi]")

ax[2].plot(tnb, spin1nb, alpha=0.66, label=r"$\Omega_1$")
ax[2].plot(tsec, spin1sec, alpha=0.66, ls="--", c="tab:blue")
ax[2].plot(tnb, spin2nb, alpha=0.66, label=r"$\Omega_2$")
ax[2].plot(tsec, spin2sec, alpha=0.66, ls="--", c="tab:orange")
ax[2].plot(tsec, n, alpha=0.66, c="black", ls=":", label="$n$")
ax[2].plot(tsec, necc, alpha=0.66, c="green", ls=":", label=r"$n_{\rm ecc}$")
ax[2].plot(tsec, nsync, alpha=0.66, c="red", ls=":", label=r"$n_{\rm p\hbox{-}s}$")
ax[2].set_ylabel("angular frequency [rad/yr]")
ax[2].set_yscale("log")
ax[2].legend(loc="lower right", ncols=2)

ax[3].plot(tnb, obl1nb, alpha=0.66, label=r"$\theta_1$")
ax[3].plot(tsec, obl1sec, alpha=0.66, ls="--", c="tab:blue")
ax[3].plot(tnb, obl2nb, alpha=0.66, label=r"$\theta_2$")
ax[3].plot(tsec, obl2sec, alpha=0.66, ls="--", c="tab:orange")
ax[3].set_ylabel("obliquity [deg]")
ax[3].set_xlabel("time / initial period")
ax[3].legend(loc="upper right")

for axx in ax:
    axx.margins(x=5e-3)
    axx.grid()

plt.show()