#! /usr/bin/env python
"""Ploting data for nucleosome crosslinking.\

Read a pickled list of interaction of 5' position in nucleosome \
base pairs against LYS residues and plot the results.\
There are more routines than needed (the ones that build who is who, copied \
from the script that generates the pickle file in the first place)
"""

import argparse
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import scipy.fftpack
import sys


ROLLO = "\n \
Process polymerase stop data (average) and compute log2 ratio\
"
whowhen = 'Guillem Portella, v0.01, 03-2017'


def smooth(x, window_len):
    """Do smooth things up."""
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    w = np.ones(window_len, 'd')
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[len(x[window_len - 1:0:-1]):len(x[window_len - 1:0:-1]) +
             len(x) + 1]


def compute_enrichment(dt, sg_window=19):
    """Compute average of two replicates, ration and log2."""
    dt["av_fc"] = dt[["cl_1", "cl_2"]].mean(axis=1)
    dt["av_ctrl"] = dt[["free_1", "free_2"]].mean(axis=1)
    # dt["log2fc"] = smooth(np.log2(dt["av_fc"] / dt["av_ctrl"]), window)
    if sg_window <= 8:
        print "Can not use this coeff polynomial order for SavGol"
        print "You get this error because window is equal or lower than 9"
        sys.exit(0)
    dt["log2fc"] = savgol_filter(
        np.log2(dt["av_fc"] / dt["av_ctrl"]), sg_window, 8, mode='mirror')
    dt["fft"] = 2.0 / len(dt["log2fc"]) * \
        np.abs(scipy.fftpack.fft(dt["log2fc"]))


def parse_arguments():
    """Parse the arguments.

    rtype: arguments dictionary
    """
    parser = argparse.ArgumentParser(description=ROLLO, epilog=whowhen)
    parser.add_argument(
        '-fwd', metavar='fwd', nargs='?', required=True,
        help='Data from the fwd strand ')
    parser.add_argument(
        '-rev', metavar='rev', nargs='?', required=True,
        help='Data from the rev strand ')
    parser.add_argument(
        '-o',
        metavar='output_f',
        nargs='?', required=False, default="poly_stop",
        help='Outputfile basename for log2 calculation')
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    """Average and proccess data from polymerase stop sequencing.

    Right now only does averaging of fC, control and computes log2(ratio)
    in fwd and reverse strand.
    The output could be fed to read_bp_list_and_plot.py as experimental data
    """
    args = parse_arguments()
    dt_fwd = pd.read_table(args["fwd"])
    dt_rev = pd.read_table(args["rev"])
    enrichment = []
    compute_enrichment(dt_fwd, sg_window=9)
    out_fwd = args["o"] + "_fwd.dat"
    out_fwd_ft = args["o"] + "_fwd_ft.dat"
    out_rev = args["o"] + "_rev.dat"
    out_rev_ft = args["o"] + "_rev_ft.dat"
    dt_fwd[["pos", "log2fc"]].to_csv(
        out_fwd, sep=" ", header=False, index=False)
    dt_fwd[["fft"]].to_csv(
        out_fwd_ft, sep=" ", header=False, index=False)
    compute_enrichment(dt_rev, sg_window=9)
    dt_rev[["pos", "log2fc"]].to_csv(
        out_rev, sep=" ", header=False, index=False)
    dt_rev[["fft"]].to_csv(
        out_rev_ft, sep=" ", header=False, index=False)
