import pandas as pd
import numpy as np
import argparse
import os

def convert_keplerian(txtfile):
    out = np.loadtxt(txtfile)
    inds = np.array(out[:, 0], dtype=np.int64)

    elems = ["a", "e", "i", "ome", "Ome", "nu"]
    cols = ["k", "Nexcurs", "Nexchan", "coll"]
    cols = cols + [x + "_inn_f" for x in elems] + [x + "_out_f" for x in elems] + [x + "_sin_f" for x in elems]
    cols = cols + ["begin_time", "end_time", "t2enc", "ibin", "jbin", "ksin"]
    cols = cols + ["i_c", "j_c", "m1_i", "m2_i", "m3_i", "Mbh"]
    cols = cols + [x + "_inn_i" for x in elems] + [x + "_out_i" for x in elems] + [x + "_sin_i" for x in elems]
    cols = cols + ["v_rel"]

    outdf = pd.DataFrame(out[:, 1:], index=inds, columns=cols)
    outdf.sort_index(inplace=True)

    h5name = txtfile.replace(".txt", ".h5")
    outdf.to_hdf(h5name, key="output")


def convert_isol(txtfile):
    out = np.loadtxt(txtfile)
    inds = np.array(out[:, 0], dtype=np.int64)

    elems = ["a", "e", "i", "ome", "Ome", "nu"]
    cols = ["k", "Nexcurs", "Nexchan", "coll"]
    cols = cols + [x + "_inn_f" for x in elems] + [x + "_out_f" for x in elems]
    cols = cols + ["begin_time", "end_time", "t2enc", "ibin", "jbin", "ksin"]
    cols = cols + ["i_c", "j_c", "m1_i", "m2_i", "m3_i"]
    cols = cols + [x + "_inn_i" for x in elems] + [x + "_out_i" for x in elems]

    outdf = pd.DataFrame(out[:, 1:], index=inds, columns=cols)
    outdf.sort_index(inplace=True)

    h5name = txtfile.replace(".txt", ".h5")
    outdf.to_hdf(h5name, key="output")

    exmapfile = txtfile.replace("_output", "_exmap")
    if os.path.exists(exmapfile):
        exmap = np.loadtxt(exmapfile)
        cols = ["fname", "nex", "tapo", "mbin", "msin", "ainn", "einn", "iinn", "Omeinn", "aout", "eout", "iout", "Omeout"]
        outdf = pd.DataFrame(exmap, columns=cols)
        outdf.set_index(['fname', 'nex'], inplace=True)
        outdf.sort_index(inplace=True)
        h5name = exmapfile.replace(".txt", ".h5")
        outdf.to_hdf(h5name, key="output")


def option_parser():
    result = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    result.add_argument("txtname", help="txt name", type=str, nargs='+')
    result.add_argument("-isol", help="isolation sim", action="store_true", dest="isol")
    return result


if __name__ == "__main__":
    args = option_parser().parse_args()
    for txtf in args.txtname:
        if args.isol:
            convert_isol(txtf)
        else:
            convert_keplerian(txtf)
