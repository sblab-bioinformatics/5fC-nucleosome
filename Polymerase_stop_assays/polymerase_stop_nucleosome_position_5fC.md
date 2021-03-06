# Position of nucleosome phase based on polymerase stop experiments for 5fC

Polymerase stop experiments upon 5fC-nucleosome ligation can give information on the actual positioning of the 5fC around the histone core. The cross-linking takes place between the formyl of the 5fC and the lysines.

The idea is to compute a probability map of sorts, indicating which bases along the Widom sequence would be close to a lysine, aiding the identification of cross-linking lysines.

## Base-pair--Lys interaction likelyhood

From a set of molecular dynamics trajectories of dinucleosome system, we first split the trajectory into individual nucleosomes, concatenate these individual trajectories, and compute the probability to find potential 5fC-lys contacts.

### Prepare the MD trajectories

The all-atoms simulations are akin to those published in Collepardo et al. JACS 2015 paper. The manuscript containts all the relevant information to replicate these simulations. 

We concatenated all the dinucleosome simulations into one large trajector, removing ions and water. The trajectories had been previously imaged. We used `-skip 10` to reduce further the number of frames.

To reduce the overhead of writing a code to deal with dinucleosomes, and since the sequences are exactly the same, we will split the trajectory into two nucleosomes and then concatenate them into one.


```bash
#! /usr/bin/env bash

source ~/GMX_2016_2/bin/GMXRC.bash

traj_folder=`pwd`/'concat_traj'
pdb_file=${traj_folder}/'init_novsites_dinuc.pdb'
xtc_file=${traj_folder}/'all_concat_novistes_dinuc.xtc'


echo 'chain A or chain B or chain C or chain D or chain E or chain F or chain G or chain H or chain I or chain J
chain K or chain L or chain M or chain N or chain O or chain P or chain Q or chain R or chain S or chain T
name 13 chains_A_to_J
name 14 chains_K_to_T
q' > kk

gmx  make_ndx -f ${pdb_file}  -o  index_split.ndx < kk && rm kk

echo 13 | gmx editconf -f ${pdb_file} -o first_nuc.pdb -n index_split.ndx -resnr 1

echo 13 | gmx trjconv -f ${xtc_file} -o first.xtc -n index_split.ndx
echo 14 | gmx trjconv -f ${xtc_file} -o second.xtc -n index_split.ndx

# Concatenate and fix time step (cosmetic)
gmx trjcat -f first.xtc second.xtc -cat -o concat.xtc &&
gmx trjconv -f concat.xtc -tu ns -timestep 2 -o concatenated_dinuc.xtc &&
rm -fr concat.xtc first.xtc second.xtc

echo ' chain A or chain B
q ' > k
echo 13 | gmx make_ndx  -f first_nuc.pdb -o index_fitting.ndx < k && rm k
echo 13 0  |  gmx trjconv -f concatenated_dinuc.xtc -n index_fitting.ndx -s first_nuc.pdb -o concat_dinuc_fit.xtc -fit rot+trans &&
rm concatenated_dinuc.k first_nuc.pdb second_nuc.pdb second.pdb index_split.ndx index_fitting.ndx
rm -fr \#*
```

### Finding contact distances

We wrote a python code,[`build_bp_lys_contacts.py`](scripts/build_bp_lys_contacts.py), you'll also need [this](scripts/nucleic_protein_db.py)), to compute the closest lysines to each possible nucleotide in a nucleosome. Check the header of the `.py` files for dependencies. From the MD trajectory above, we centered the distance calculation to the position that the 5-formyl group would adopt, given that that residue was a C. Then we keep a list of dictionaries for each base pair position (1 to 147 for both strands). Each element in the list is a dictionary that stores the time-averaged contact with those lysines that are less than a cutoff (set to 1.2nm). The contacts are made continuous via a switching function (saturation-like curve) such that anything below 0.5nm gives a value of 1, and anything above that decays to 0 (at ~1 nm is almost zero).

```python
%matplotlib inline
import matplotlib.pyplot as plt
import math
import numpy as np

def switch2(b, d0, x):
    return 1.0 / (1 + math.exp(b*(x - 1.5*d0)))

p = np.arange(0, 2, 0.01)
s = np.zeros(200)
for i in range(200):
    s[i] = switch2(10, 0.5, float(i)/100)

fig = plt.figure()
fig.suptitle('Switching function 1.0 / (1.0 + exp(b*(x - 1.5*d0))) ', fontsize=14)
plt.title("b=10, d0=0.5")
plt.ylabel("f(d)", fontsize=20)
plt.xlabel("d", fontsize=20)
plt.plot(p, s)
plt.savefig("switch_function.png")
```

![switch_function](pictures/switch_function.png)

The keys of the dictionary are the residue number of each lysine, and the values are the time-averaged contacts. The code also computes the phase profile along the nucleosome, to get and indicator of which groove is facing the core. The way the phase is defined now is the angle (in rad) between a vector pointing outwards from the center of mass of a base pair \[actually the N1 of the pyr] towards the minor groove. So a small "phase" value associated with a given base pair indicates that the minor groove is facing the histone core.

A second code ([`read_bp_list_and_plot.py`](scripts/read_bp_list_and_plot.py)) takes care to plot the data contained as a pickled file written by `build_bp_lys_contacts.py`.


With `./XXX.py -h` you get the options. The source code is relatively well documented. A typical executions would go like this,

```bash
#! /usr/bin/env bash

# one with the tails; takes long
./build_bp_lys_contacts.py -f concat_dinuc_fit.xtc -s first_nuc_res.pdb -st 4252 -ol list_bp_lys_dict_tails.pkl

# one without the tails
./build_bp_lys_contacts.py -f concat_dinuc_fit.xtc -s first_nuc_res.pdb -st 4252 -notails -ol list_bp_lys_dict_notails.pkl
```

## Analysis including new data set (fwd and reverse strands)

We now have both strands, plus controls. The control in case of 5fC crosslinked against the histones is the experiment done without the crosslinking.

Using [do_process_fwd_rev_data.py](scripts/do_process_fwd_rev_data.py).

```bash

exe_f=~/work/repos/sblab_projects/projects/20151125_5fC_nucleosome/20170130_polymerase_stop_modelling/scripts/do_process_fwd_rev_data.py
${exe_f} -fwd 20170315_table_forward.txt -rev 20170315_table_reverse.txt -o poly_stop

# add the data to the plot

exe_f=read_bp_list_and_plot.py
#no tails
${exe_f} -s first_nuc_res.pdb -lys list_bp_lys_dict_notails.pkl  -phase bp_phase_profile.txt -exp poly_stop_fwd.dat poly_stop_rev.dat -o lys_contact_map_notails
# tails
${exe_f} -s first_nuc_res.pdb -lys list_bp_lys_dict_tails.pkl  -phase bp_phase_profile.txt -o lys_contact_map_tails -exp poly_stop_fwd.dat poly_stop_rev.dat
```

We used this [notebook](scripts/plot_polstop_profile.ipynb) to plot the data. 

# Overlap analysis between polstop data and simulations

Used [read_bp_list_process_likelihood.py](scripts/read_bp_list_process_likelihood.py), e.g.,

```bash
exe_f=read_bp_list_process_likelihood.py
${exe_f} -s first_nuc_res.pdb -lys list_bp_lys_dict_tails.pkl  -phase bp_phase_profile.txt -exp poly_stop_fwd.dat poly_stop_rev.dat  -seq widom601.fa
```

After analysis of the data, here a table linking nucleotide positions with lysine residues. The phase indicates minor groove (negative) and major groove (positive).

Combining modelling high-probabilty regions with polstop high-stop regions (product)

| DNA res. | Near K res. |  Phase [-1,1] | Histone | Hist. tail? |
| -------- | -------------- |  ----- | -------| ------- |
| 27 | 661 (0.37) | -0.56 | H2b_1 | Yes |
| 41 | 506 (0.35) | 0.71 | H4_1 |  -  |
| 46 | 1158 (0.18) | -0.83 | H2b_2 | Yes |
| 61 | 460 (0.62) | 0.37 | H4_1 |  -  |
| 61 | 460 (0.62) | 0.37 | H4_1 |  -  |
| 71 | 896 (0.30) | 0.29 | H3_2 |  -  |
| 79 | 790 (0.09) | -0.79 | H3_2 | Yes |
| 91 | 845 (0.45) | 0.03 | H3_2 |  -  |
| 170 (124) | 1170 (0.26) | 0.24 | H2b_2 | Yes |
| 174 (120) | 1031 (0.27) | -0.56 | H2a_2 | Yes |
| 188 (106) | 993 (0.31) | 0.71 | H4_2 |  -  |
| 193 (101) | 684 (0.29) | -0.83 | H2b_1 | Yes |
| 208 (86) | 947 (0.60) | 0.37 | H4_2 |  -  |
| 208 (86) | 947 (0.60) | 0.37 | H4_2 |  -  |
| 218 (76) | 409 (0.31) | 0.29 | H3_1 |  -  |
| 226 (68) | 659 (0.04) | -0.79 | H2a_1 |  -  |
| 230 (64) | 460 (0.12) | 0.90 | H4_1 |  -  |
| 238 (56) | 358 (0.41) | 0.03 | H3_1 |  -  |
| 242 (52) | 1154 (0.08) | 0.69 | H2b_2 | Yes |
| 242 (52) | 1154 (0.08) | 0.69 | H2b_2 | Yes |
| 253 (41) | 1167 (0.04) | 0.39 | H2b_2 | Yes |

Using only polymerase stop data (all of them, not only the most significant 
in terms of logFC(5fC/control).

| DNA res. | Near K res. |  Phase [-1,1] | Histone | Hist. tail? |
| -------- | -------------- |  ----- | -------| ------- |
| 22 | 683 (0.02) | 0.64 | H2b_1 | Yes |
| 27 | 661 (0.37) | -0.56 | H2b_1 | Yes |
| 35 | 1154 (0.04) | -0.37 | H2b_2 | Yes |
| 35 | 1154 (0.04) | -0.37 | H2b_2 | Yes |
| 41 | 506 (0.35) | 0.71 | H4_1 |  -  |
| 44 | 1154 (0.05) | 0.18 | H2b_2 | Yes |
| 50 | 1155 (0.03) | 0.36 | H2b_2 | Yes |
| 53 | 449 (0.02) | 0.75 | H4_1 |  -  |
| 61 | 460 (0.62) | 0.37 | H4_1 |  -  |
| 65 | 437 (0.02) | 0.21 | H3_1 |  -  |
| 69 | 1137 (0.12) | -0.51 | H2a_2 |  -  |
| 79 | 790 (0.09) | -0.79 | H3_2 | Yes |
| 83 | 947 (0.11) | 0.90 | H4_2 |  -  |
| 91 | 845 (0.45) | 0.03 | H3_2 |  -  |
| 95 | 932 (0.04) | 0.69 | H4_2 | Yes |
| 106 | 667 (0.45) | 0.39 | H2b_1 | Yes |
| 106 | 667 (0.45) | 0.39 | H2b_1 | Yes |
| 115 | 1031 (0.06) | 0.84 | H2a_2 | Yes |
| 118 | 1031 (0.44) | -0.31 | H2a_2 | Yes |
| 124 | 1170 (0.31) | 0.63 | H2b_2 | Yes |
| 139 | 308 (0.25) | -0.32 | H3_1 | Yes |
| 139 | 308 (0.25) | -0.32 | H3_1 | Yes |
| 156 (138) | 350 (0.22) | 0.60 | H3_1 |  -  |
| 156 (138) | 350 (0.22) | 0.60 | H3_1 |  -  |
| 158 (136) | 1092 (0.20) | 0.80 | H2a_2 |  -  |
| 164 (130) | 1154 (0.14) | -0.48 | H2b_2 | Yes |
| 170 (124) | 1170 (0.26) | 0.24 | H2b_2 | Yes |
| 174 (120) | 1031 (0.27) | -0.56 | H2a_2 | Yes |
| 182 (112) | 1023 (0.03) | -0.37 | H2a_2 | Yes |
| 182 (112) | 1023 (0.03) | -0.37 | H2a_2 | Yes |
| 188 (106) | 993 (0.31) | 0.71 | H4_2 |  -  |
| 191 (103) | 921 (0.07) | 0.18 | H4_2 | Yes |
| 193 (101) | 684 (0.29) | -0.83 | H2b_1 | Yes |
| 199 (95) | 932 (0.09) | 0.99 | H4_2 | Yes |
| 202 (92) | 924 (0.00) | 0.03 | H4_2 | Yes |
| 208 (86) | 947 (0.60) | 0.37 | H4_2 |  -  |
| 212 (82) | 785 (0.06) | 0.21 | H3_2 | Yes |
| 216 (78) | 649 (0.14) | -0.51 | H2a_1 |  -  |
| 220 (74) | 409 (0.55) | 0.95 | H3_1 |  -  |
| 228 (66) | 460 (0.28) | 0.22 | H4_1 |  -  |
| 230 (64) | 460 (0.12) | 0.90 | H4_1 |  -  |
| 238 (56) | 358 (0.41) | 0.03 | H3_1 |  -  |
| 242 (52) | 1154 (0.08) | 0.69 | H2b_2 | Yes |
| 253 (41) | 1167 (0.04) | 0.39 | H2b_2 | Yes |
| 253 (41) | 1167 (0.04) | 0.39 | H2b_2 | Yes |
| 262 (32) | 661 (0.03) | 0.84 | H2b_1 | Yes |

