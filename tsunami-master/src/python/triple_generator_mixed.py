import numpy as np
from tsunami import stateid
from triple_generator_common import okinami_option_parser, process_args
from triple_generator_nbody import TripleSystemNbody
from triple_generator_secular import TripleSystemSecular


class TripleSystemMixed():
    id = stateid

    def __init__(self,
                 m1, m2, m3, a1, a2, e1, e2, i_mut, ome1, ome2, Ome1, R1, R2, R3,
                 k1=0.0, k2=0.0, tau1=0.0, tau2=0.0,
                 nocol=False, tide=False, gr=True, filename="test_triple", mute_output=False,
                 stabcrit=1,
                 nopn1=False,
                 nopn2=False,
                 nopn25=False,
                 check_quasi_secular=True):

        self.filename = filename
        self.filename_secular = filename + "_sec"
        self.filename_nbody = filename + "_nb"
        self.TripleSec = TripleSystemSecular(m1=m1, m2=m2, m3=m3, a1=a1, a2=a2, e1=e1, e2=e2, i_mut=i_mut,
                                             ome1=ome1, ome2=ome2, Ome1=Ome1, R1=R1, R2=R2, R3=R3,
                                             k1=k1, k2=k2, tau1=tau1, tau2=tau2, nocol=nocol, tide=tide, gr=gr,
                                             filename=self.filename_secular, mute_output=mute_output,
                                             nopn1=nopn1, nopn2=nopn2, nopn25=nopn25, stabcrit=stabcrit,
                                             check_quasi_secular=check_quasi_secular)

    def run(self, ftime, dtout, logtime=False, timing=True):
        self.TripleSec.run(ftime, dtout, logtime=logtime, timing=timing)
        print("Exited secular run")
        if self.TripleSec.secularcode.hascollided:
            print("Run has collided, closing")
            return

        if self.TripleSec.time >= ftime:
            print("Run has finished, closing")
            return

        if self.TripleSec.killed:
            print("Run was killed, closing")
            return

        if not self.TripleSec.secularcode.isdynstable:
            print("Run became unstable, continuing")

        if self.TripleSec.secularcode.stopsec:
            print("Run became quasi-secular, continuing")

        #self.TripleSec.plot(task=args.task)
        e1, e2, ome1, ome2, Ome1, a1, H, a2, m1, m2, m3, R1, R2, R3 = self.TripleSec.y
        k1, k2 = self.TripleSec.k1, self.TripleSec.k2
        tau1, tau2 = self.TripleSec.tau1, self.TripleSec.tau2
        nu1 = np.random.uniform(0.0, 2*np.pi)
        nu2 = np.random.uniform(0.0, 2*np.pi)
        begin_time = self.TripleSec.time

        self.TripleNb = TripleSystemNbody(m1=m1, m2=m2, m3=m3, a1=a1, a2=a2, e1=e1, e2=e2, i_mut=i_mut,
                                          ome1=ome1, ome2=ome2, Ome1=Ome1, R1=R1, R2=R2, R3=R3,
                                          k1=k1, k2=k2, tau1=tau1, tau2=tau2, nu1=nu1, nu2=nu2,
                                          filename=self.filename_nbody, begin_time=begin_time)
        self.TripleNb.run(ftime, dtout, timing=timing)



if __name__ == "__main__":
    args = okinami_option_parser().parse_args()

    m1, m2, m3, R1, R2, R3, a1, e1, a2, e2, i_mut, ome1, ome2, Ome1, k1, k2, tau1, tau2 = process_args(args)

    kzt = TripleSystemMixed(m1, m2, m3,
                              a1, a2,
                              e1, e2,
                              i_mut,
                              ome1, ome2,
                              Ome1,
                              R1, R2, R3,
                              k1, k2,
                              tau1, tau2,
                              args.nocol, args.tide, args.gr,
                              nopn25=args.nopn25, nopn2=args.nopn2, nopn1=args.nopn1,
                              check_quasi_secular=args.qsec,
                              stabcrit=args.stabcrit,
                              filename=args.filename,
                              mute_output=args.mute_output)

    kzt.run(args.ftime, args.dt_out)