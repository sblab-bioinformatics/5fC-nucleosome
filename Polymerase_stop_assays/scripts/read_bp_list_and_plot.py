#! /usr/bin/env python
"""Ploting data for nucleosome crosslinking.\

Read a pickled list of interaction of 5' position in nucleosome \
base pairs against LYS residues and plot the results.\
There are more routines than needed (the ones that build who is who, copied \
from the script that generates the pickle file in the first place)
"""


import mdtraj as md
import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import nucleic_protein_db as nap_db
import pickle
import datetime
from Bio import SeqIO
import bz2
from itertools import islice
import peakutils

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

histones = {"H3_1": [295, 429], "H4_1": [439, 531],
            "H2a_1": [532, 659], "H2b_1": [660, 781],
            "H3_2": [782, 916], "H4_2": [917, 1018],
            "H2a_2": [1019, 1146], "H2b_2": [1147, 1168]}

histones_tails = {"H3_1": [295, 329], "H4_1": [439, 444],
                  "H2a_1": [532, 539], "H2b_1": [660, 746],
                  "H3_2": [782, 817], "H4_2": [917, 939],
                  "H2a_2": [1019, 1034], "H2b_2": [1147, 1175]}


class basepair():
    """Document.
    """

    def __init__(self):
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
        """Check because I think I compute more things than needed.
        """
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


def window(seq, n=2):
    "Return a sliding window (of width n) over data from the iterable."
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    """obviously from http://stackoverflow.com/a/6822773"""
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
    """Builds an array of basepair objects
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
    """Explain."""
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
    """Builds two dictionaries and a list
    the first one lys residue number as key, and a list of sidechains as
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


def rescale_zero_to_one(phase):
    max_v = phase[:, 1].max()
    min_v = phase[:, 1].min()
    d = max_v - min_v
    return (phase[:, 1] - min_v) / d


def read_fasta(args):
    """Read the sequences."""
    if args['seq'].split(".")[-1] == "bz2":
        handle = bz2.BZ2File(args['f'], "r")
    else:
        handle = open(args['seq'], "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    return records


def plot_data(list_dict, phase, experiment, sequence, peaks,
              outpng="plot.png"):
    """Do plot the data."""
    title = "Contacts maps"
    plt.title(title)
    num_plots = 1
    n_bp = len(phase)
    contact = np.zeros((2, n_bp))
    x_ax = np.arange(1, n_bp + 1)
    f, ax = plt.subplots(num_plots, 1, sharex=True,
                         figsize=(18, 4), dpi=120)
    ax2 = ax.twiny()
    f.subplots_adjust(hspace=0.2)
    ax.grid(b=True, which='major',
            color='gray', linestyle='--', linewidth=0.5)
    plt.locator_params(nbins=20)
    ax.set_axisbelow(True)
    strand = 0
    [i.set_linewidth(0.5) for i in ax.spines.itervalues()]
    for i, d in enumerate(list_dict):
        # we append on a different array if we pass n_bp limit
        if i >= n_bp:
            j = i - n_bp
            strand = 0
        else:
            j = i
            strand = 1
        max_v = find_closest(d)
        contact[strand, j] = max_v
    ax.plot(x_ax, phase,
            label="Phase", color="black", linestyle="-", linewidth=0.5)
    ax.plot(x_ax, contact[0, :],
            label="First strand", color="red", linestyle="-")
    ax.plot(x_ax, contact[1, :],
            label="Second strand", color="blue", linestyle="-")
    ax.plot(experiment[0][:, 0], experiment[0][:, 1],
            label="Polym. stop fwd", color="orange", linestyle="--")
    ax.plot(experiment[1][:, 0], experiment[1][:, 1],
            label="Polym. stop rev", color="blueviolet", linestyle="--")

    ax.set_ylabel("Signal (for now)")
    ax.set_xlabel("Base pair number")
    ax.legend(ncol=4).get_frame().set_linewidth(0.1)

    # print the sequence on top
    ax2.set_xlim(ax.get_xlim())
    ax2_tick_loc = np.arange(NUC_LEN)
    ax2.set_xticks(ax2_tick_loc)
    ax2.set_xticklabels([b for b in sequence])
    ax2.set_xlabel(r"Sequence fwd strand")
    # color all Cs
    # [i.set_color("blue") if i.get_text() ==
    # 'C' else '' for i in ax2.get_xticklabels()]
    # color only CGs
    for ii in window(ax2.get_xticklabels()):
        if ii[0].get_text() == "C" and ii[1].get_text() == "G":
            ii[0].set_color("red"), ii[1].set_color("red")

    plt.savefig(outpng)


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
        '-op',
        metavar='lys_contact_map',
        nargs='?', required=False, default="lys_contact_map.png",
        help='Distance to closest lysines as function of base postion')
    args = vars(parser.parse_args())
    return args


def read_data(pick_file, prof_file, exp):
    """Read all the data and return the numpy arrays.

    The phase data gets rescaled from 0 to 1, for the sake of it
    """
    with open(pick_file, "r") as lh:
        list_dict = pickle.load(lh)
    phase_prof = np.genfromtxt(prof_file)
    rescaled_phase = rescale_zero_to_one(phase_prof)
    experiment = []
    experiment.append(np.genfromtxt(exp[0]))
    experiment.append(np.genfromtxt(exp[1]))
    return list_dict, rescaled_phase, experiment


def do_call_peaks(polstop):
    """Uses peakutils module to find the peak location."""
    indexes = []
    for exp in polstop:
        ind = peakutils.indexes(exp[:, 1], thres=0.2, min_dist=2)
        indexes.append(ind)
    return indexes


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

    # outpng = "_".join([str(datetime.date.today().strftime("%Y%m%d")),
    #                   arguments["op"]])
    outpng = arguments["op"]
    lys_dict, phase_prof, experiment = read_data(ldict, phase_data, exp)
    peaks = do_call_peaks(experiment)
    plot_data(lys_dict, phase_prof, experiment, sequence, peaks, outpng)
