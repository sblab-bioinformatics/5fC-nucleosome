#! /usr/bin/env python3
"""Submits deeptools compute matrix."""

import os
import shutil
import tarfile
import slurm
from contextlib import contextmanager
from collections import OrderedDict
import sys


def make_tarfile(output_filename, source_dir):
    """Make a tarfile."""
    with tarfile.open(output_filename, "w:bz2") as tar:
        tar.add(source_dir, arcname=os.path.basename(source_dir))
        return


def make_dir(dir_name):
    """Make a new dir if not there."""
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    else:
        print("File was already there")


def remove_dir(dir_name):
    """Remove directory if there."""
    if os.path.exists(dir_name):
        shutil.rmtree(dir_name)


@contextmanager
def cd(newdir):
    """Safely cd, from http://stackoverflow.com/a/24176022 ."""
    prevdir = os.getcwd()
    print(prevdir)
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

class referencepoint(object):
    """Same as below, but center."""

    def __init__(self, *initial_data, **kwargs):
        """Init."""
        self.len_after = 1000
        self.reference = "Reference"
        self.score_list = ["score"]
        self.out_matrix = "matrix.mat.gz"
        self.len_before = 1000
        self.len_after = 1000
        # the input line
        self.inp_line = """
        computeMatrix reference-point \
        -R {reference} \
        -S {score_list} \
        -out {out_matrix} \
        --referencePoint "center"\
        -b {len_before} \
        -a {len_after} \
        -bs 10 \
        --skipZeros \
        --sortRegions "no" \
        -p "max" """

        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        """Write output if print() is called on object."""
        options = self.__dict__
        return self.inp_line.format(**options)

    def cmd(self):
        """Work like __str__."""
        options = self.__dict__
        return self.inp_line.format(**options)

class centerplotprofile(object):
    """Prepare input for plotprofile.

    Inherits methods from scaleregions
    """

    def __init__(self, *initial_data, **kwargs):
        """Init."""
        self.matrix_file = "matrix.mat.gz"
        self.out_plot = "plot.pdf"
        self.out_dat = "plot.dat"
        self.plot_type = "se"
        self.colors = ["red"]
        self.refpoint_label = "center"
        self.samples_label = ["A1"]
        self.regions_label = "Region"
        self.y_label = "RPKM"
        # the input line
        self.inp_line = """
        plotProfile \
        --matrixFile {matrix_file} \
        -out {out_plot} \
        --outFileNameData {out_dat} \
        --plotType "{plot_type}" \
        --colors {colors} \
        --samplesLabel {samples_label} \
        --regionsLabel {regions_label} \
        -y "{y_label}" \
        --perGroup  """

        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        """Write output if print() is called on object."""
        options = self.__dict__
        return self.inp_line.format(**options)

    def cmd(self):
        """Work like __str__."""
        options = self.__dict__
        return self.inp_line.format(**options)


class plotprofile(object):
    """Prepare input for plotprofile.

    Inherits methods from scaleregions
    """

    def __init__(self, *initial_data, **kwargs):
        """Init."""
        self.matrix_file = "matrix.mat.gz"
        self.out_plot = "plot.pdf"
        self.out_dat = "plot.dat"
        self.plot_type = "se"
        self.colors = ["red"]
        self.start_label = "start"
        self.end_label = "end"
        self.samples_label = ["A1"]
        self.regions_label = "Region"
        self.y_label = "RPKM"
        # the input line
        self.inp_line = """
        plotProfile \
        --matrixFile {matrix_file} \
        -out {out_plot} \
        --outFileNameData {out_dat} \
        --plotType "{plot_type}" \
        --colors {colors} \
        --startLabel "{start_label}" \
        --endLabel "{end_label}" \
        --samplesLabel {samples_label} \
        --regionsLabel {regions_label} \
        -y "{y_label}" \
        --perGroup  """

        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        """Write output if print() is called on object."""
        options = self.__dict__
        return self.inp_line.format(**options)

    def cmd(self):
        """Work like __str__."""
        options = self.__dict__
        return self.inp_line.format(**options)


class scaleregions(object):
    """Keep track of all the deeptools options.

    The members can be initialized from a dict. If ther
    are not there the member will be created.
    What we should actually do is look for specific
    entries in the dictionary, if they are not there then
    use a default, and ignore the weird entries.
    This should probably have been a function, but I
    need to refresh my OO abilities.
    """

    def __init__(self, *initial_data, **kwargs):
        """Init."""
        self.len_after = 1000
        self.reference = "Reference"
        self.score_list = ["score"]
        self.out_matrix = "matrix.mat.gz"
        self.start_label = "start"
        self.start_label = "end"
        self.region_len = 500
        self.len_before = 1000
        self.len_after = 1000
        # the input line
        self.inp_line = """
        computeMatrix scale-regions \
        -R {reference} \
        -S {score_list} \
        -out {out_matrix} \
        -m {region_len} \
        --startLabel "{start_label}" \
        --endLabel "{end_label}" \
        -b {len_before} \
        -a {len_after} \
        -bs 10 \
        --skipZeros \
        --sortRegions "no" \
        -p "max" """

        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        """Write output if print() is called on object."""
        options = self.__dict__
        return self.inp_line.format(**options)

    def cmd(self):
        """Work like __str__."""
        options = self.__dict__
        return self.inp_line.format(**options)


# Some setup
path_bw = "/scratcha/sblab/martin03/repository/"\
    + "20151125_5fC_nucleosome/data/20170209/bam_clean"
path_bed_5fC = "/scratcha/sblab/portel01/"\
    + "project/euni_5fc/intersect_histone_marks_cpgi/"

marks = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3", "H3K27me3", "CpGi"]

marks_bed_KO_name = path_bed_5fC + "{:s}_consensus_KO.bed.gz"
marks_bed_WT_name = path_bed_5fC + "{:s}_consensus_WT.bed.gz"
marks_bed_WT_KO_name = path_bed_5fC + "{:s}_common_WT_KO.bed.gz"

# dictionaries where key is a mark name and the values is an array with
# the name of the _5fC (first entry), and the _no5fC (second) file name
bed_marks_WT = dict((mark, [marks_bed_WT_name.format(mark + "_5fC"),
                            marks_bed_WT_name.format(mark + "_no5fC")]) for mark in marks)
bed_marks_KO = dict((mark, [marks_bed_KO_name.format(mark + "_5fC"),
                            marks_bed_KO_name.format(mark + "_no5fC")]) for mark in marks)
bed_marks_WT_KO = dict((mark, [marks_bed_WT_KO_name.format(mark + "_5fC"),
                               marks_bed_WT_KO_name.format(mark + "_no5fC")]) for mark in marks)


# The MNase data
WT_A1 = 'ear030_HindbrainA1.150618.ncbi37.bw'
WT_C1 = 'ear031_HindbrainC1.150618.ncbi37.bw'
WT_F1 = 'ear032_HindbrainF1.150618.ncbi37.bw'
KO_K1 = 'ear040_K1_hindbrain.SLX-10219.ncbi37.bw'
KO_K2 = 'ear041_K2_hindbrain.SLX-10219.ncbi37.bw'
bw_files = OrderedDict((name, path_bw + "/" + eval(name))
                       for name in ['WT_A1', 'WT_C1', 'WT_F1', 'KO_K1', 'KO_K2'])

# FC pulldown --> regions where 5fC were found
KO_2A2 = 'grm064_TDG_KO_2A2_hindbrain.clean.bw'
KO_O4 = 'grm065_TDG_KO_O4_hindbrain.clean.bw'
KO_Q2 = 'grm066_TDG_KO_Q2_hindbrain.clean.bw'
KO_Z1 = 'grm067_TDG_KO_Z1_hindbrain.clean.bw'
WT_2A8 = 'grm068_TDG_WT_2A8_hindbrain.clean.bw'
WT_O3 = 'grm069_TDG_WT_O3_hindbrain.clean.bw'
WT_Q3 = 'grm070_TDG_WT_Q3_hindbrain.clean.bw'
WT_Z3 = 'grm071_TDG_WT_Z3_hindbrain.clean.bw'
invivo_bw_files = OrderedDict((name, eval(name))
                              for name in ['WT_A1', 'WT_C1',
                                           'WT_F1', 'KO_K1', 'KO_K2'])


conditions = {"WT": bed_marks_WT, "KO": bed_marks_KO, "WT_KO": bed_marks_WT_KO}
# Start mark iteration
for condition, dict_condition in conditions.items():
    for mark, mark_bed_file in dict_condition.items():

        # using templated printing options
        # mark_bed_file should allow for a list as well

        options_center_reg = {"reference": " ".join(mark_bed_file),
                             "score_list": " ".join(bw_files.values()),
                             "out_matrix": mark + "_" + condition +  "_referencepoint.mat.gz",
                             "region_len": 500,
                             "len_before": 1000,
                             "len_after": 1000
                             }

        options_centerplotprofile = {"matrix_file": mark + "_" + condition + "_referencepoint.mat.gz",
                               "out_plot": mark + "_" + condition +  "_referencepoint.pdf",
                               "out_dat": mark + "_" + condition + "_referencepoint.dat",
                               "plot_type": "se",
                               "colors": " ".join(["red", "red", "red", "blue", "blue"]),
                               "samples_label": " ".join(bw_files.keys()),
                               "regions_label": " ".join(['"' + mark + '"',
                                                          '"' + mark + '_no5fC"']),
                               "y_label": "RPKM"
                               }

        refpoint_reg = referencepoint(options_center_reg)
        plot_reg = centerplotprofile(options_centerplotprofile)
        make_dir(condition+"_center")
        with cd(condition+"_center"):
            make_dir(mark)
            with cd(mark):
                with open("submit_" + mark + ".sh", "w") as out:
                    cmd = refpoint_reg.cmd() + " && \ " + plot_reg.cmd()
                    out.write(cmd)
                    slurm.sbatch(cmd, J=mark, log=mark + ".log")
