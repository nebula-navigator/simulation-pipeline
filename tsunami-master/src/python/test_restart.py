import tsunami
import numpy as np
import matplotlib.pyplot as plt

def run(restart=False, dt=0.1, ft=65, restart_time=25):
    code = tsunami.Tsunami(1, 1)
    m = np.array([3., 4., 5., 6.])
    p = np.array([[1.,3.,0.], [-2.,-1.,0.], [1.,-1.,0.], [2.,0.,0.]])
    v = np.array([[0.,0.,0.], [0.,0.,0.], [0.,0.,0.], [0.,0.,0.]])
    r = np.array([0., 0., 0., 0.])
    st = np.array([0, 1, 2, 3])
    code.add_particle_set(p, v, m, r, st)
    code.sync_internal_state(p, v)
    totp = [p.copy()]

    # Obtain chain vectors (only a copy, not the real thing)
    chp, chv = code.get_chain_vectors()
    print(chp, chv)

    restarted = False
    realt = 0
    data = []

    while realt < ft:
        realt = realt + dt

        code.evolve_system(realt)
        realt = code.time

        if restart and not restarted and realt > restart_time:
            print("Stopping code and restarting at time", realt)

            chp, chv = code.get_chain_vectors()
            timestep = code.timestep
            
            masses = np.arange(4, dtype=np.float64)
            code.sync_masses(masses)
            print("Masses", masses)

            print("Chain pos:\n", chp)
            print("Chain vel:\n", chv)
            print("timestep", timestep)
            print("dtphysical", code.dtphysical)

            code.save_restart_file("pytest_ics.bin")
            code.load_restart_file("pytest_ics.bin")

            chp_f, chv_f = code.get_chain_vectors()
            timestep_f = code.timestep

            print("Chain pos:\n", chp_f)
            print("Chain vel:\n", chv_f)
            print("timestep", timestep_f)
            print("dtphysical", code.dtphysical)

            print("Are chain positions the same:", np.all(chp_f == chp))
            print("Are chain velocities the same:", np.all(chv_f == chv))
            print("Is timestep the same:", np.all(timestep == timestep_f))
            del chp_f, chv_f

            restarted = True

        # Synchronize coordinates to interface
        code.sync_internal_state(p, v)

        data.append([code.time, code.deltaE, code.energy])
        totp.append(p.copy())
        print("time={:g}/{:g} - {:2.1f}% - dE/E0={:e} - E={:g}".format(realt, ft, realt/ft*100,
                                                                       code.deltaE, code.energy), end="\033[K\r")

        #if restarted: break

    totp = np.vstack(totp)
    data = np.asarray(data)

    return data, totp

rstime = 25
data_nores, pos_nores = run(restart=False, restart_time=rstime)
data_res, pos_res= run(restart=True, restart_time=rstime)

f, ax = plt.subplots(2, figsize=(6,6), sharex="all", gridspec_kw=dict(hspace=0),  tight_layout=True)

ax[0].plot(data_nores[:,0], data_nores[:,1], lw=1.5)
ax[1].plot(data_nores[:,0], data_nores[:,2], lw=1.5)
ax[0].plot(data_res[:,0], data_res[:,1], lw=1.5)
ax[1].plot(data_res[:,0], data_res[:,2], lw=1.5)
ax[0].set_ylabel("$\Delta E / E_0$")
ax[1].set_ylabel("$E$")
ax[1].set_xlabel("time")
for axx in ax: axx.axvline(rstime, c='grey', ls="--")

plt.show()


f, ax = plt.subplots(1, figsize=(6,6),  tight_layout=True)

Nt = len(pos_nores)
for i in range(4):
    ax.plot(pos_nores[i:Nt-i:4,0], pos_nores[i:Nt-i:4,1], lw=1.5)
    ax.plot(pos_res[i:Nt-i:4,0], pos_res[i:Nt-i:4,1], lw=1.5, ls=":")

plt.show()
