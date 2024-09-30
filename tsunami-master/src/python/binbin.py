import tsunami
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib as mpl

mpl.rcParams['figure.raise_window'] = False

keplutils = tsunami.KeplerUtils()

a1 = a2 = 1 # au
e1 = e2 = 0.0
m1 = m2 = m3 = m4 = 10 # MSun
a3 = -3.2 # au
e3 = 2

print("Peri:", a3*(1.0-e3), "au")

nu1, nu2 = np.pi, 0.6*np.pi

d3 = 1000 # au
maxth = np.arccos(-1 / e3)
latum = a3 * (1 - e3 * e3)
nu3 = latum / d3 - 1
nu3 = nu3 / e3
nu3 = - np.arccos(nu3)

print("Initial dist:", d3, "au")

i1 = i2 = i3 = 0.0
ome1 = ome2 = ome3 = 0.0
Ome1 = Ome2 = Ome3 = 0.0

pos_vel2 = np.array([0.,0.,0., 0.,0.,0.])
minn1 = m1+m2
pos_vel1 = keplutils.kepl_to_cart(pos_vel2, m1, m2, a1, e1, i1, ome1, Ome1, nu1)
pos_vel_inn_com = (pos_vel1*m1 + pos_vel2*m2)/minn1
pos_vel1 -= pos_vel_inn_com
pos_vel2 -= pos_vel_inn_com

pos_vel4 = np.array([0.,0.,0., 0.,0.,0.])
minn2 = m3+m4
pos_vel3 = keplutils.kepl_to_cart(pos_vel4, m3, m4, a2, e2, i2, ome2, Ome2, nu2)
pos_vel_inn_com = (pos_vel3*m3 + pos_vel4*m4)/minn2
pos_vel3 -= pos_vel_inn_com
pos_vel4 -= pos_vel_inn_com

mout = m1+m2+m3+m4
pos_vel_b1 = np.array([0.,0.,0., 0.,0.,0.])
pos_vel_b2 = keplutils.kepl_to_cart(pos_vel_b1, minn1, minn2, a3, e3, i3, ome3, Ome3, nu3)

pos_vel1 += pos_vel_b1
pos_vel2 += pos_vel_b1
pos_vel3 += pos_vel_b2
pos_vel4 += pos_vel_b2

quadruple_com = (pos_vel1*m1 + pos_vel2*m2 + pos_vel3*m3 + pos_vel4*m4)/mout
pos_vel1 -= quadruple_com
pos_vel2 -= quadruple_com
pos_vel3 -= quadruple_com
pos_vel4 -= quadruple_com

code = tsunami.Tsunami()

m = np.array([m1, m2, m3, m4])
print("m", m)

p = np.array([pos_vel1[:3], pos_vel2[:3], pos_vel3[:3], pos_vel4[:3]])
print("p", p)

v = np.array([pos_vel1[3:], pos_vel2[3:], pos_vel3[3:], pos_vel4[3:]])
print("v", v)

r = np.array([0., 0., 0., 0.])
st = np.array([-1, -1, -1, -1])
code.add_particle_set(p, v, m, r, st)
N = code.N

a1_l, e1_l, t_l = [], [], []
a2_l, e2_l = [], []
a3_l, e3_l = [], []

hfont = {'fontname':'DejaVu Serif', 'fontsize' : 19, 'labelpad' : 10 }

fig = plt.figure(figsize=(10, 7))

plt.ioff()

gs = gridspec.GridSpec(2, 2)
ax1a = plt.subplot(gs[0, 0])
ax1b = plt.subplot(gs[1, 0])
ax2 = plt.subplot(gs[:,1])
ax2.set_aspect('equal')
ax2.tick_params(labelsize=14)

#plt.tight_layout()

ax1a.set_ylabel("e", **hfont)
ax1a.tick_params(labelsize=14)
ax1b.set_ylabel("i [deg]", **hfont)
ax1b.set_xlabel("time [yr]", **hfont)
ax1b.tick_params(labelsize=14)
#f.subplots_adjust(hspace=0)

totp = p

time = 0.0
ftime = 5e2 * 2*np.pi
fdist = d3
dtout = 2*np.pi
while time < ftime:
    time = time + dtout
    
    code.evolve_system(time)
    time = code.time
    code.sync_internal_state(p, v)
    t_l.append(time/(2*np.pi))
    totp = np.vstack((totp, p))

    ax2.cla()
    ax1a.cla()
    ax1b.cla()

    #for ip in range(N):
    ax2.set_xlim(-70, 70)
    ax2.set_ylim(-d3/9, d3/9)
    ax2.scatter(p[:,0], p[:,1], s=50, color=['tab:blue', 'tab:orange', 'tab:green', 'tab:red'])

    for ip in range(N):
        ax2.plot(totp[ip::N,0], totp[ip::N,1], lw=5.5)
          #print(ax2.get_xlim(), ax2.get_ylim())

    orb1 = keplutils.cart_to_kepl(p[0]-p[1], v[0]-v[1], m1, m2)
    p_inn1 = (p[0]*m1 + p[1]*m2)/minn1
    v_inn1 = (v[0]*m1 + v[1]*m2)/minn1

    orb2 = keplutils.cart_to_kepl(p[2]-p[3], v[2]-v[3], m3, m4)
    p_inn2 = (p[2]*m3 + p[3]*m4)/minn2
    v_inn2 = (v[2]*m3 + v[3]*m4)/minn2

    orb3 = keplutils.cart_to_kepl(p_inn1-p_inn2, v_inn1-v_inn2, minn1, minn2)

    a1_l.append(orb1[0])
    a2_l.append(orb2[0])
    a3_l.append(orb3[0])
    e1_l.append(orb1[1])
    e2_l.append(orb2[1])
    e3_l.append(orb3[1])

    ax1a.plot(t_l, a1_l, lw=2.5)
    ax1a.plot(t_l, a2_l, lw=2.5)
    ax1b.plot(t_l, e1_l, label='binary 1', lw=2.5)
    ax1b.plot(t_l, e2_l, label='binary 2', lw=2.5)
    ax1b.legend()
    ax1a.set_ylabel("semimajor axis [au]")
    ax1b.set_xlabel("time [yr]")
    ax1b.set_ylabel("eccentricity")

    ax1b.set_ylim(0, 1)
    ax1a.set_ylim(0.0, max(a1_l[0]+1, max(a1_l)))

    plt.draw()
    plt.pause(1e-3)

    dist = p_inn1-p_inn2
    dist = (dist*dist).sum()**0.5

    if dist > fdist:
        print("Reached final distance\033[K")
        break

    print("Time (direct) = {:3.2f} yr".format(time/(2*np.pi)), end="\r")

plt.show()
