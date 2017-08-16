#!/usr/bin/env python
__doc__ = """
Functionality: compute time averaged contacts between 5 position in nucleosome base pairs
against LYS residues. Position 5 as if it was a pyrimidine.
Reads: a nuclosme trajectory in XTC format and a corresponding PDB file
Outputs: a list of dictionaries in python pickle format, a phase profile for each
base pair position. The list of dictionaries contains an entry for each base position
(2*number of base pairs), and each dictionary entry has a LYS residue index as key
and its time averaged contacts as values.
To compute contacts we use a cutoff of 1.2, and a switching function that more or
less gives 1 if the distance between any sides chain atom of LYS is closer to 0.5, and
decreases to 0 like a sigmoid curve.
"""


import numba
import mdtraj as md
import sys
import argparse
import numpy as np
import math
import progressbar
import nucleic_protein_db as nap_db
import pickle


ROLLO = "\n \
Probabilities of nucleosome-lys interaction.\
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
    """Store information on the base pair.

    Basically the index (atom number) of atoms in a base
    pair. A class is probably too much, as we actually
    don't have methods.
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
        # self.c_w = 0


@numba.jit(nopython=True)
def bp_direction(bpc, wn, cn, xyz):
    """Compute the orientation of the base pair minor groove
    as the angle bisecting the vectors from the ~COM of the base pair
    and atoms on the Watson and Crick strand. The way it is defined
    here we get a vector pointing outwards from the minor groove.
    """
    vect_c_w = xyz[0, bpc, :] - xyz[0, wn, :]
    vect_c_c = xyz[0, bpc, :] - xyz[0, cn, :]
    mod_c_w = np.sqrt((vect_c_w * vect_c_w).sum())
    mod_c_c = np.sqrt((vect_c_c * vect_c_c).sum())
    sum_v = vect_c_w / mod_c_w + vect_c_c / mod_c_c
    mod_sum_v = np.sqrt((sum_v * sum_v).sum())
    return sum_v / mod_sum_v


def inrange(val, range_list):
    """Check if value is bewteen two value, both included."""
    if val >= range_list[0] and val <= range_list[1]:
        return True
    else:
        return False


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


def histcore_com(x, frame_ind=0):
    """Compute the center of mass of the histone core.

    rtype: float, the center of mass.
    """
    com = np.zeros(3)
    count = 0
    for a in x.topology.atoms:
        if a.name not in nap_db.nucleic_type_dict:
            com = com + x.xyz[frame_ind, a.index, :]
            count = count + 1
    return com / count


def return_N1_index(x_top, r_ind_1, r_ind_2):
    """For a given residue pair, return the N1 of the purine.
    rtype: int, the atom index.
    """
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
        bps.append(bp)

    return bps


def profile_phase(com_his, bps, x):
    """ Compute the angle between the vector pointing from the COM
    of the histones to the 'center' of the base pair, and the vector
    pointing from the 'center' of the base pair towards the minor groove.
    This is done for each base pair. With this definition, the angle is
    maximum when the major groove points towards the histone core.
    rtype: numpy array of phase angles.
    """
    phase = np.zeros(len(bps))
    for i, bp in enumerate(bps):
        bp_dir = bp_direction(bp.bpc, bp.wn, bp.cn, x.xyz)
        cc = x.xyz[0, bp.bpc, :] - com_his
        mod = np.sqrt((cc * cc).sum(axis=0))
        cc = cc / mod
        # both vectors are already normalized
        phase[i] = np.arccos(np.clip(np.dot(cc, bp_dir), -1.0, 1.0))
    return phase


def build_lys_list(pdb_x, b_notails):
    """Builds two dictionaries and a list
    the first one lys residue number as key, and a list of sidechains as
    value. The second has lys residue number as key and the chain number
    as value. The list is are indices all non - histone tail lys sidechain
    atoms.
    """
    lys_chain_dict = {}
    lys_at_dict = {}
    all_lys = []
    if b_notails:
        for i, c in enumerate(pdb_x.topology.chains):
            for a in c.atoms:
                if a.residue.name == "LYS" and a.name in lys_chain\
                    and len(np.where(np.searchsorted(hist_tail_boundaries,
                                                     a.residue.index,
                                                     side="right") % 2)[0]) == 0:
                    lys_chain_dict.setdefault(a.residue.index, i)
                    lys_at_dict.setdefault(a.residue.index, []).append(a.index)
                    all_lys.append(a.index)
    else:
        for i, c in enumerate(pdb_x.topology.chains):
            for a in c.atoms:
                if a.residue.name == "LYS" and a.name in lys_chain:
                    lys_chain_dict.setdefault(a.residue.index, i)
                    lys_at_dict.setdefault(a.residue.index, []).append(a.index)
                    all_lys.append(a.index)

    return [lys_chain_dict, lys_at_dict, all_lys]


def get_referencepos_b(x, res_num):
    """Compute the reference position for the 5 position
    If the base is pur, get the reference position to be the COM
    of the N7 with O6 or N6.
    If the base is pyr, get either the H5 or the C7
    rtype: np array dimension 3 (xyz coordinates)
    """
    ref_position = np.zeros(3)
    if nap_db.nucleic_type_dict[x.topology.residue(res_num).name] == "pu":
        for a in x.topology.residue(res_num).atoms:
            # this should find the COM between N7 and O6|N6
            if a.name == "N7" or a.name == "O6" or a.name == "N6":
                ref_position += x.xyz[0, a.index, :]
        ref_position = ref_position / 2
    else:
        for a in x.topology.residue(res_num).atoms:
            # this should find the COM between N7 and O6|N6
            if a.name == "H5" or a.name == "C7":
                ref_position = x.xyz[0, a.index, :]

    return ref_position


@numba.jit(nopython=True)
def compute_dis2(x, val, ref):
    """Compute the distance squared between an atom in
    x.xyz and a reference 3D vector
    The function is trivial and could be just place directly
    in the find_lys_distance code, but this way I can easily
    decorate it with numba.jit(nopython=True) and get a 4x speed up.
    rtype: float, a distance squared
    """
    return ((x[0, val, :] - ref) *
            (x[0, val, :] - ref)).sum()


def find_lys_distance(x, cutoff, lys_info, ref_position, counts):
    """Iterates over all lysines in lys_info[1], which is the
    the dictionary of side chain atoms of all lysines.

    For each lysine
    side chain determines the minimum distance within a cutoff wrt to a
    reference position (index of the W or C atom). The distance is then
    converted to a contact from [0,1] (1 being a close contact) using
    a switching function.
    """
    sq_cutoff = cutoff * cutoff
    # would be nice to make it a set
    for ind in lys_info[1]:
        min_v = 1000
        found = False
        for val in lys_info[1][ind]:
            d2 = compute_dis2(x, val, ref_position)
            if sq_cutoff > d2:
                d = math.sqrt(d2)
                found = True
                if d < min_v:
                    min_v = d
        if found:
            counts[ind] = counts.setdefault(ind, 0) + switch_f(min_v)
        found = False
        min_v = 1000


def compute_contacts(x, bp_array, lys_info, mega):
    """Compute the contacts between bases 5pyr position and lys residues.

    Updates the list of dictionaries 'mega'.
    """
    cutoff = 1.2
    n_bp = len(bp_array)
    for i, bp in enumerate(bp_array):
        # we do the W strand
        ref_position_w = get_referencepos_b(x, bp.w)
        find_lys_distance(x.xyz, cutoff, lys_info, ref_position_w, mega[i])
        j = 2 * n_bp - i - 1
        # now the C strand
        ref_position_c = get_referencepos_b(x, bp.c)
        find_lys_distance(x.xyz, cutoff, lys_info, ref_position_c, mega[j])


def switch_f(x, b=10, d0=0.5, l0=1.5):
    """Switching function

    The function is 1 between 0 and d0, and decays
    to zero at around 1.
    """
    return 1.0 / (1 + math.exp(b * (x - l0 * d0)))


def crazy_f(x, a=6, b=1, d0=0.1, r0=1):
    """Explain."""
    return (1 + (2**(a / b) - 1) * ((x - d0) / r0)**a)**(-b / a)


def init_c(bp_n):
    """Creat a list of dictionaries.
    The list containts as many entries as bases in the nucleosome.
    """
    list_d = []
    for bp in xrange(2 * bp_n):
        d = {}
        list_d.append(d)
    return list_d


def parse_arguments():
    """Parse the arguments.
    rtype: arguments dictionary
    """
    parser = argparse.ArgumentParser(description=ROLLO, epilog=whowhen)
    parser.add_argument(
        '-s', metavar='pdb', nargs='+', required=True, help='A pdb file ')
    parser.add_argument(
        '-f',
        metavar='xtc',
        nargs='+',
        required=True,
        help='An xtc file to go with the pdb ')
    parser.add_argument(
        '-n',
        metavar='n_bp',
        type=int,
        nargs='?',
        default='147',
        help='Number of DNA base pairs')
    parser.add_argument(
        '-ol',
        metavar='pickled_list_of_dict',
        nargs='?', required=False, default="list_bp_lys_dict.pkl",
        help='A pickled list containing for each base pair a dictionary of \
    lysines residues and their distance to the base C5-like position')
    parser.add_argument(
        '-op',
        metavar='phase_profile',
        nargs='?', required=False, default="bp_phase_profile.txt",
        help='The profile of the bp phase along the nucleosome')
    parser.add_argument(
        '-st',
        metavar='steps_xtc',
        type=int,
        nargs='?',
        default='1000',
        help='Number of steps to read from xtc')
    parser.add_argument(
        '-notails',
        default=False,
        action='store_true',
        help='Do not take the histone tails into account')
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    arguments = parse_arguments()
    n_bp = arguments["n"]
    last_step = arguments["st"]
    b_notails = arguments["notails"]
    bar = progressbar.ProgressBar(max_value=last_step)
    pdb_x = md.load(arguments["s"])
    lys_info = build_lys_list(pdb_x, b_notails)
    bp_array = build_base_pairs(n_bp, pdb_x)
    prof = np.zeros(arguments["n"])
    # initialize list of dictionaries
    mega = init_c(arguments["n"])
    cutoff = 1.2
    for steps, chunk in enumerate(md.iterload(arguments["f"][0],
                                              top=arguments["s"][0], chunk=1)):
        com_hist = histcore_com(chunk, 0)
        prof += profile_phase(com_hist, bp_array, chunk)
        compute_contacts(chunk, bp_array, lys_info, mega)
        bar.update(steps + 1)
        if steps >= last_step - 1:
            break

    # Write the dicts out AND normalize for storage
    for base_pos, m in enumerate(mega):
        if len(m) > 0:
            print base_pos + 1, " : ",
            for k in m:
                m[k] = m[k] / (cutoff * (steps + 1))
                print k, m[k],
            print
    # Write the profile
    with open(arguments["op"], "w") as profh:
        for i, p in enumerate(prof):
            print >> profh,  i + 1, p / (steps + 1)
    # Write the list of dictionaries to a pickled file
    picklename = arguments["ol"]
    with open(picklename, 'wb') as ph:
        pickle.dump(mega, ph)
