import argparse
import os
import numpy as np
import subprocess
import re


class SubmitTsunamiJob():
    """
    This is the script I use to submit long simulations on HPC clusters using the PBS queue system.
    It expects 1 single run for 1 single job on 1 single core.
    If you want to launch multiple short simulations on many single core jobs, use the ChunkTsunamiJob class
    It expect that the set of simulations has the following structure:
    ./setname_folder/
                    0000001.dat
                    0000002.dat
                    .......
                    "{:07d}".format(Nmax)
    where each 7-digits file is Tsunami initial condition
    """

    ic_re = re.compile("(\d+).dat")

    def __init__(self,
                 ftime,
                 dtime,
                 walltime,
                 folder,
                 ijob=0,
                 fjob=-1,
                 Npart=None,
                 wTides=None,
                 wPNs=None,
                 ptype=None,
                 oflags="",
                 Nmax=False,
                 test=False):

        folder = os.path.normpath(folder)
        self.folder = folder
        self.setname = os.path.basename(folder)
        self.wh = int(walltime)
        self.wm = round(60 * (walltime % 1))
        self.ftime = ftime  # *(2 * np.pi)*1e6
        self.dtime = dtime  # *(2 * np.pi)*1e6

        self.test = test

        # Flags
        self.N = Npart
        self.wTides = wTides
        self.wPNs = wPNs
        self.ptype = ptype
        self.oflags = oflags
        self.cont = True
        self.moduleload = ["gcc/8.2.0"]  # add here modules to be loaded
        self.username = "tranias"  # add here your user name
        self.Nmax = Nmax

        self.icfiles = [x for x in os.listdir(self.folder) if self.ic_re.match(x)]
        self.icfiles.sort()
        self.Nic = len(self.icfiles)
        print("Found {:d} icfiles".format(self.Nic))
        self.icfiles = self.icfiles[ijob:fjob]
        print("Selected {:d} (i->j {:d}->{:d})".format(len(self.icfiles), ijob, fjob))

        # Checks which jobs are running, assuming they were launched with the submitter. Works with PBS
        with open("joblist.dsv", 'w') as jbf:
            subprocess.run(["qstat", "-u", self.username, "-a", "-w", "-f", "-F", "dsv"], stdout=jbf)

        jobtab = np.loadtxt("joblist.dsv", delimiter='|', usecols=(0, 1), dtype=('S100', 'S100'))

        self.joblis = [x.decode().replace("Job_Name=", "") for x in jobtab[:, 1]] if len(jobtab) != 0 else []

    def add_module(self, module):
        self.moduleload.append(module)

    def save_options(self):
        pass

    def load_options(self):
        pass

    def check_and_submit(self):
        self.Nrest, self.Nbegin, self.Nskip = 0, 0, 0
        for icf in self.icfiles:
            colf = icf.replace(".dat", "_collision.dat")
            colf = os.path.join(self.folder, colf)
            if os.path.isfile(colf):
                continue
            if self.Nmax is not None and self.Nrest + self.Nbegin >= self.Nmax:
                print("Reached {:d} submitted jobs, quitting".format(self.Nmax))
                print("  (stopped at {:s})".format(icf))
                break
            enf = icf.replace(".dat", "_energy.dat")
            enf = os.path.join(self.folder, enf)
            if os.path.isfile(enf):
                with open(enf, "r") as enef:
                    lines = enef.readlines()
                    lastene = np.loadtxt(lines[-1:])
                    self.lastt = lastene[0]
                if self.lastt < self.ftime:
                    if self.is_running(icf): continue
                    self.make_pbs_job(icf)
                    self.Nrest += 1
            else:
                self.lastt = 0
                if self.is_running(icf): continue
                self.make_pbs_job(icf)
                self.Nbegin += 1

        print("Nrest", self.Nrest)
        print("Nbegin", self.Nbegin)
        print("Nskip", self.Nskip)

    def is_running(self, icf):
        inds = icf.replace(".dat", "")
        self.jobname = self.setname + "_n{:d}".format(int(inds))
        if self.jobname in self.joblis:
            self.Nskip += 1
            if self.test: subprocess.run(["echo", "skipping " + self.jobname + ", is running"])
            return True
        return False

    def make_pbs_job(self, icf):
        inds = icf.replace(".dat", "")
        tmpfile = "tmp_" + self.setname + "_" + inds + ".pbs"
        launchline = "./tsunami.x {:s}/{:s}.dat -ft {:1.3e} -dt {:1.3e}".format(self.folder, inds, self.ftime,
                                                                                self.dtime)

        tmpjob = open(tmpfile, "w")
        tmpjob.write("#!/bin/bash\n")
        tmpjob.write("#PBS -N {:s}\n".format(self.jobname))
        tmpjob.write("#PBS -q mid\n")
        tmpjob.write("#PBS -l walltime={:d}:{:02d}:00\n".format(self.wh, self.wm))
        tmpjob.write("\n")
        tmpjob.write("cd $PBS_O_WORKDIR\n")
        for mod in self.moduleload:
            tmpjob.write("module load {:s}".format(mod))
        tmpjob.write("\n")

        if self.N:
            launchline = launchline + " -N {:d}".format(self.N)
        if self.ptype:
            pts = ["{:d}".format(x) for x in self.ptype]
            launchline = launchline + " -pt " + " ".join(pts)
        if self.wPNs:
            launchline = launchline + " -PN"
        if self.wTides:
            launchline = launchline + " -T"
        if self.cont:
            launchline = launchline + " -c"
        if self.oflags:
            launchline = launchline + self.oflags

        tmpjob.write(launchline)

        tmpjob.write("\n")
        tmpjob.close()

        if (not self.test):
            subprocess.run(["qsub", tmpfile])
        else:
            subprocess.run(["echo", "qsub " + tmpfile + "  t={:1.3e} ".format(self.lastt) + launchline])

        os.remove(tmpfile)


def new_argument_parser():
    result = argparse.ArgumentParser()
    # exclude = result.add_mutually_exclusive_group()
    result.add_argument("folder",
                        help="folder", type=str)
    result.add_argument("-i", dest="ijob", default=0,
                        type=int, help="in code units")
    result.add_argument("-f", dest="fjob", default=-1,
                        type=int, help="in code units")
    result.add_argument("-w", dest="walltime", default=6,
                        type=float, help="in hours")

    result.add_argument("-ft", dest="ftime", default=2 * np.pi * 1e6,
                        type=float, help="in code units")
    result.add_argument("-dt", dest="dtime", default=2 * np.pi * 1e2,
                        type=float, help="in code units")
    result.add_argument("-T", dest="wTides", action="store_true",
                        help="tides", default=None)
    result.add_argument("-PN", dest="wPNs", action="store_true",
                        help="postnewtonians", default=None)
    result.add_argument("-pt", dest="ptype", nargs="*", default=None,
                        type=int, help="particle type")
    result.add_argument("-N", dest="Npart", default=None,
                        type=int, help="N particles")
    result.add_argument("-F", dest="oflags", default="",
                        type=str, help="other flags")

    result.add_argument("-Nmax", dest="Nmax", default=None,
                        help="max number of jobs", type=int)
    result.add_argument("--test", dest="test", action="store_true",
                        help="test")
    return result


if __name__ in ('__main__'):
    args = vars(new_argument_parser().parse_args())
    PBS = SubmitTsunamiJob(**args)

    # PBS.make_pbs_job("0000000")
    PBS.check_and_submit()
