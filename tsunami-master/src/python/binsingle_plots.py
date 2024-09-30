import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py


def traj(fname, task):
    # All is in N-body units (G=1, L=au, M=MSun, T=yr/2pi)
    h5f = h5py.File(fname, 'r')
    pos = h5f["p"][()]
    # vel = h5f["v"][()]
    # time = h5f["t"][()]
    h5f.close()

    f, ax = plt.subplots(1, 1, tight_layout=True)
    N = 3
    Nt = len(pos) // N

    ic, jc = 0, 1

    for ib in range(N):
        s = np.s_[ib:Nt * N - N + ib:N]
        ax.plot(pos[s, ic], pos[s, jc], alpha=0.66)
        ax.scatter(pos[s, ic][0], pos[s, jc][0], marker="o", alpha=0.66)

    if "show" in task:
        plt.show()


def run_options():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result.add_argument("fname", help="filename", type=str)
    result.add_argument("-k", help="type of plot",
                        dest="kplot", type=str, default="traj")
    result.add_argument("-t", "-task", help="save format. options: show, png, pdf",
                        dest="task", type=str, nargs="+", default=["show"])
    return result


if "__main__" == __name__:
    args = run_options().parse_args()

    if args.kplot == "traj":
        traj(args.fname, args.task)
