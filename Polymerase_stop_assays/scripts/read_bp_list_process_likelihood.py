#! /usr/bin/env python
"""Combine polstop and simulations data for nucleosome crosslinking.

Read a pickled list of interaction of 5' position in nucleosome
base pairs against LYS residues and the experimental data
and output the most likey interaction partners in positions of
significant polymerase stop.
There are more routines than needed (the ones that build who is who, copied
from the script that generates the pickle file in the first place)
"""


import mdtraj as md
import math
import sys
import argparse
import numpy as np
import nucleic_protein_db as nap_db
import pickle
from Bio import SeqIO
import bz2
from itertools import islice
import peakutils
import scipy.optimize
import operator

NUC_LEN = 147

ROLLO = "\n \
Read the nucleosome-lys interaction list of dictionaries\
"
whowhen = 'Guillem Portella, v0.01, 02-2017'

# HG1, HD1, etc... get renamed to HG3, HD3, etc... in mdtraj
lys_chain = ["CG", "HG3", "HG2", "CD", "HD3", "HD2",
             "CE", "HE3", "HE2", "NZ", "HZ1", "HZ2", "HZ3"]

# conversion from histo_tails to a list containing [,) pairs, where
# even indeces are entry points, and odd indeces are exit point
# http://stackoverflow.com/questions/23639361/fast-checking-of-ranges-in-python
hist_tail_boundaries = [295, 330, 439, 445, 523, 540,
                        660, 747, 782, 818, 917, 940, 1019, 1035, 1147, 1176]
hist_boundaries = [295, 439, 439, 532, 532, 660,
                   660, 782, 782, 917, 917, 1019, 1019, 1146, 1147, 1269]
# TOBEDONE
index_hist_boundaries = {0: "H3_1", 1: "H4_1", 2: "H2a_1", 3: "H2b_1",
                         4: "H3_2", 5: "H4_2", 6: "H2a_2", 7: "H2b_2"}

histones = {"H3_1": [295, 438], "H4_1": [439, 531],
            "H2a_1": [532, 659], "H2b_1": [660, 781],
            "H3_2": [782, 916], "H4_2": [917, 1018],
            "H2a_2": [1019, 1146], "H2b_2": [1147, 1268]}

histones_tails = {"H3_1": [295, 329], "H4_1": [439, 444],
                  "H2a_1": [532, 539], "H2b_1": [660, 746],
                  "H3_2": [782, 817], "H4_2": [917, 939],
                  "H2a_2": [1019, 1034], "H2b_2": [1147, 1175]}


class basepair():
    """Document."""

    def __init__(self):
        """Initialize."""
        # each has a watson and crick strand
        # the id is the residue number
        self.w = 0
        self.c = 0
        # each strand has its reference N1/N9, index is the atom index
        # plus there we set the center of the bp as the N1 of pyr
        # this should work fine for what we need, but it can be sophisticated
        self.wn = 0  # N1/N9 of watson
        self.cn = 0  # N1/N9 of crick
        self.bpc = 0  # we set the center of the bp as the N1 of pyr
        self.c_w = 0

    def bp_normal(self, x, frame_ind=0):
        """Check because I think I compute more things than needed."""
        vect_c_w = x.xyz[frame_ind, self.bpc, :] - x.xyz[frame_ind, self.wn, :]
        vect_c_c = x.xyz[frame_ind, self.bpc, :] - x.xyz[frame_ind, self.cn, :]
        # norm_vect = np.cross(vect_c_w, vect_c_c)
        mod_c_w = np.sqrt((vect_c_w * vect_c_w).sum(axis=0))
        mod_c_c = np.sqrt((vect_c_c * vect_c_c).sum(axis=0))
        sum_v = vect_c_w / mod_c_w + vect_c_c / mod_c_c
        mod_sum_v = np.sqrt((sum_v * sum_v).sum(axis=0))
        return sum_v / mod_sum_v


def inrange(val, range_list):
    """Document todo."""
    if val >= range_list[0] and val <= range_list[1]:
        return True
    else:
        return False


def fit_sin(tt, yy):
    """Fit sin to the input time sequence, and return fitting parameters.

    "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"
    from http://stackoverflow.com/a/42322656
    Originally it returns many things, but i am only interested in the fit
    """
    tt = np.array(tt)
    yy = np.array(yy)
    # assume uniform spacing
    ff = np.fft.fftfreq(len(tt), (tt[1] - tt[0]))
    Fyy = abs(np.fft.fft(yy))
    # excluding the zero frequency
    guess_freq = abs(ff[np.argmax(Fyy[1:]) + 1])
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array(
        [guess_amp, 2. * np.pi * guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c): return A * np.sin(w * t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
#    f = w / (2. * np.pi)

    # def fitfunc(t): return A * np.sin(w * t + p) + c
    return A * np.sin(w * tt + p) + c
    # return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f,
    #       "period": 1. / f, "fitfunc": fitfunc, "maxcov": np.max(pcov),
    #       "rawres": (guess, popt, pcov)}


def window(seq, n=2):
    """Move a windo-w (of width n) over data from the iterable.

    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    obviously from http://stackoverflow.com/a/6822773
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def return_N1N9_index(x_top, res_ind):
    """Return index of N1/N9 atoms.

    For a given residue index, returns the atom index of the N9
    or N1 nitrogen, if it's a NA
    """
    if nap_db.nucleic_type_dict[x_top.residue(res_ind).name] == "pu":
        for x in x_top.residue(res_ind).atoms:
            if x.name == "N9":
                atom_index = x.index
    elif nap_db.nucleic_type_dict[x_top.residue(res_ind).name] == "py":
        for x in x_top.residue(res_ind).atoms:
            if x.name == "N1":
                atom_index = x.index
    else:
        print "Error, residue name not in the database"
        print x_top.residue(res_ind).name
        sys.exit(0)

    return atom_index


def return_N1_index(x_top, r_ind_1, r_ind_2):
    """For a given residue pair, return the N1 of the purine."""
    if nap_db.nucleic_type_dict[x_top.residue(r_ind_1).name] == "pu":
        for x in x_top.residue(r_ind_1).atoms:
            if x.name == "N1":
                atom_index = x.index
    if nap_db.nucleic_type_dict[x_top.residue(r_ind_2).name] == "pu":
        for x in x_top.residue(r_ind_2).atoms:
            if x.name == "N1":
                atom_index = x.index
    return atom_index


def build_base_pairs(n_bp, pdb_x):
    """Build an array of basepair objects.

    rtype: a list of basepair objects
    """
    x_top = pdb_x.topology
    bps = []
    for i in xrange(n_bp):
        j = 2 * n_bp - i - 1
        bp = basepair()
        bp.w = i
        bp.c = j
        bp.wn = return_N1N9_index(x_top, i)
        bp.cn = return_N1N9_index(x_top, j)
        bp.bpc = return_N1_index(x_top, i, j)
        # bp.bp_normal(pdb_x)
        bps.append(bp)

    return bps


def profile(com_his, bps, x):
    """Calculate the phase along the sequence."""
    phase = np.zeros(len(bps))
    for i, bp in enumerate(bps):
        bp_norm = bp.bp_normal(x)
        cc = x.xyz[0, bp.bpc, :] - com_his
        mod = np.sqrt((cc * cc).sum(axis=0))
        cc = cc / mod
        # both vectors are already normalized
        phase[i] = np.arccos(np.clip(np.dot(cc, bp_norm), -1.0, 1.0))
    return phase


def build_lys_list(pdb_x):
    """Build two dictionaries and a list.

    The first one lys residue number as key, and a list of sidechains as
    value. The second has lys residue number as key and the chain number
    as value. The list is are indices all non - histone tail lys sidechain
    atoms.
    """
    lys_chain_dict = {}
    lys_at_dict = {}
    all_lys = []
    for i, c in enumerate(pdb_x.topology.chains):
        for a in c.atoms:
            if a.residue.name == "LYS" and a.name in lys_chain \
                and len(np.where(np.searchsorted(hist_tail_boundaries,
                                                 a.residue.index,
                                                 side="right") % 2)[0]) == 0:
                lys_chain_dict.setdefault(a.residue.index, i)
                lys_at_dict.setdefault(a.residue.index, []).append(a.index)
                all_lys.append(a.index)

    return [lys_chain_dict, lys_at_dict, all_lys]


def find_closest(d):
    """Explain."""
    max = 0
    for x in d:
        if d[x] > max:
            max = d[x]
    return max


def find_lys_closest(d):
    """Explain."""
    max = 0
    x_max = None  # Not sure if None or False
    for x in d:
        if d[x] > max:
            max = d[x]
            x_max = x
    return x_max


def find_two_lys_closest(d):
    """Sort the dict by value, return list of tuples.

    rtype: returns a list of two tuples
    """
    sorted_tuple_d = sorted(
        d.items(), key=operator.itemgetter(1), reverse=True)
    return sorted_tuple_d[:2]


def rescale_zero_to_one(phase):
    """Do it."""
    max_v = phase[:, 1].max()
    min_v = phase[:, 1].min()
    d = max_v - min_v
    return (phase[:, 1] - min_v) / d


def rescale_min_one_to_one(phase):
    """Do it."""
    new_min = -1
    new_max = 1
    d_new = new_max - new_min
    max_v = phase[:, 1].max()
    min_v = phase[:, 1].min()
    d = max_v - min_v
    return (d_new * (phase[:, 1] - min_v)) / d + new_min


def boxify(phase):
    """Convert the phase profile to box [-1,1]."""
    return np.sign(phase[:, 1] - phase[:, 1].mean())


def read_fasta(args):
    """Read the sequences."""
    if args['seq'].split(".")[-1] == "bz2":
        handle = bz2.BZ2File(args['f'], "r")
    else:
        handle = open(args['seq'], "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    return records


def parse_arguments():
    """Parse the arguments.

    rtype: arguments dictionary
    """
    parser = argparse.ArgumentParser(description=ROLLO, epilog=whowhen)
    parser.add_argument(
        '-s', metavar='pdb', nargs='+', required=True, help='A pdb file ')
    parser.add_argument(
        '-lys',
        metavar='pklfile',
        nargs='?',
        required=True,
        help='An pickle file with lists of dictironaries, generated\
         by base_pair.py ')
    parser.add_argument(
        '-phase',
        metavar='phase_data',
        nargs='?',
        required=True,
        help='The phase of the nucleosome as function of bp')
    parser.add_argument(
        '-exp',
        nargs=2,
        metavar=("fwd_strand", "fwd_rev"),
        required=True,
        help='The data from the polymerase stop experiment, both strands')
    parser.add_argument(
        '-seq',
        metavar='seq',
        nargs='?',
        default='sequence.fa',
        help='Number of DNA base pairs')
    parser.add_argument(
        '-n',
        metavar='n_bp',
        type=int,
        nargs='?',
        default=NUC_LEN,
        help='Number of DNA base pairs')
    parser.add_argument(
        '-o',
        metavar='outfile',
        nargs='?', required=False, default="analysis_polstop.txt",
        help='Output of the program')
    args = vars(parser.parse_args())
    return args


def read_data(pick_file, prof_file, exp):
    """Read all the data and return the numpy arrays.

    The phase data gets rescaled from 0 to 1, for the sake of it
    """
    with open(pick_file, "r") as lh:
        list_dict = pickle.load(lh)
    phase_prof = np.genfromtxt(prof_file)
    # sinusoidal = fit_sin(phase_prof[:, 0], phase_prof[:, 1])
    # sign_phase = np.sign(sinusoidal - sinusoidal.mean())
    scaled_phase = rescale_min_one_to_one(phase_prof)
    experiment = []
    experiment.append(np.genfromtxt(exp[0]))
    experiment.append(np.genfromtxt(exp[1]))
    return list_dict, scaled_phase, experiment


def do_compute_contacts(lys_dict):
    """Compute a profile showing where contacts are maximum."""
    contact = np.zeros((2, n_bp))
    for i, d in enumerate(lys_dict):
        # we append on a different array if we pass n_bp limit
        if i >= n_bp:
            j = i - n_bp
            strand = 0
        else:
            j = i
            strand = 1
        max_v = find_closest(d)
        contact[strand, j] = max_v
    return contact


def do_pad_exp_and_overlap(exp, contacts):
    """Overlap (multiply) experimental and computed curves."""
    # add padding to exp to make arrays of same size
    # the fwd strans lacks data at the beginnig, and the rev at the end
    init_fwd = int(exp[0][0, 0]) - 1
    end_rev = NUC_LEN - int(exp[1][-1, 0])
    fwd_exp_padded = np.concatenate((np.zeros(init_fwd), exp[0][:, 1]))
    rev_exp_padded = np.concatenate((exp[1][:, 1], (np.zeros(end_rev))))
    overlap = []
    conv_1 = fwd_exp_padded * contacts[0, :]
    overlap.append(conv_1)
    conv_2 = rev_exp_padded * contacts[0, :]
    overlap.append(conv_2)

    padded_exp = []
    padded_exp.append(fwd_exp_padded)
    padded_exp.append(rev_exp_padded)

    return padded_exp, overlap


def do_call_peaks(x):
    """Use peakutils module to find the peak location."""
    indexes = []
    for exp in x:
        ind = peakutils.indexes(exp, thres=0.25, min_dist=1)
        indexes.append(ind)
    return indexes


def do_find_candidates(lys_dict, signal, phase):
    """Find the location of the closest C to a given signal peak.

    Reports back the lysine index of the closest C site of an signal peak.
    By signal I mean the place where polstat and/or my modelling both show a
    signal.
    """
    peaks = do_call_peaks(signal)
    results = []
    for i, peak_strand in enumerate(peaks):
        for j, p in enumerate(peak_strand):
            start = p - 3
            end = p + 4
            positions = np.array([k for k, ltr in enumerate(
                sequence[start:end]) if ltr == "C"])
            if positions.any():
                pos = positions[np.argmin(np.abs(positions - 3))]
                c_position = start + pos
                # we want a positive signal, over a little made-up threshold
                if signal[i][c_position] > 0.05:
                    # care: lys_dict was built with both strand concatenated
                    # (see do_compute_contacts)
                    pos_in_dict = c_position + i * n_bp
                    partial_res = []
                    partial_res.append(pos_in_dict)
                    who = find_two_lys_closest(lys_dict[pos_in_dict])
                    for el in who:
                        partial_res.append(el)
                    p = phase[c_position]
                    partial_res.append(p)
                    results.append(partial_res)
    return results


def print_table(results):
    """Print the results in nice markdown table."""
    # header
    print "| DNA res. | Near K res. |  Phase [-1,1] | Histone | Hist. tail? |"
    print "| -------- | -------------- |  ----- | -------| ------- | "
    for l in results:
        where = int(math.floor(np.searchsorted(hist_boundaries, l[1][0],
                                               side="right") / 2))
        who = index_hist_boundaries[where]
        if len(np.where(np.searchsorted(hist_tail_boundaries,
                                        l[1][0],
                                        side="right") % 2)[0]) != 0:
            tail = "Yes"
        else:
            tail = " - "

        if l[0] < NUC_LEN:
            print "| {:d} | {:d} ({:3.2f}) | {:3.2f} | {:s} | {:s} | ".\
                format(l[0], l[1][0], l[1][1], l[3], who, tail)
        else:
            print "| {:d} ({:d}) | {:d} ({:3.2f}) | {:3.2f} | {:s} | {:s} | ".\
                format(l[0], l[0] - NUC_LEN, l[1][0], l[1][1], l[3], who, tail)


def print_pymol(results):
    """Print the results for pymol."""
    # header
    print "select  C_candidates, ",
    for i, l in enumerate(results):
        if i != len(results) - 1:
            print "resid ", l[0] + 1, "or",
        else:
            print "resid ", l[0] + 1

    print "select  Lys_candidates, ",
    for i, l in enumerate(results):
        if i != len(results) - 1:
            print "resid ", l[1][0] + 1, "or",
        else:
            print "resid ", l[1][0] + 1


if __name__ == "__main__":
    """Plot the base-pair vs lys contacts as function of base pair position.

    Read the base-pair vs lysine contact dictionary, the phase profile and the
    experimental data. Put everything into a plot and print out as png.
    """

    arguments = parse_arguments()
    n_bp = arguments["n"]
    pdb_x = md.load(arguments["s"])
    # [lys_ch_dict, lys_at_dict, lys_at] = build_lys_list(pdb_x)
    lys_info = build_lys_list(pdb_x)
    bp_array = build_base_pairs(n_bp, pdb_x)
    ldict = arguments["lys"]
    phase_data = arguments["phase"]
    exp = arguments["exp"]
    seqrec = read_fasta(arguments)
    if len(seqrec) != 1:
        print "Please provide only one sequence"
        sys.exit(0)
    if len(seqrec[0].seq) != NUC_LEN:
        print "Hei, the length of the sequence should be ", NUC_LEN
        sys.exit(0)
    else:
        sequence = seqrec[0].seq

    lys_dict, phase_prof, experiment = read_data(ldict, phase_data, exp)
    contacts = do_compute_contacts(lys_dict)
    # we define the overalp as the element wise product of the curves
    padded_exp, overlap = do_pad_exp_and_overlap(experiment, contacts)
    results = do_find_candidates(lys_dict, overlap, phase_prof)
    print_table(results)
    print
    print "Results for pymol"
    print
    print_pymol(results)
    results_polstop = do_find_candidates(lys_dict, padded_exp, phase_prof)
    print
    print "With only polstop"
    print
    print_table(results_polstop)
    print
    print "Results for pymol"
    print
    print_pymol(results_polstop)
