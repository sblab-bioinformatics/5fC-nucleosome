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
        -bs 1 \
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
path_hist_bw = "/scratcha/sblab/portel01/project/euni_5fc/chipseq_histmarks"
path_bed_5fC = "/scratcha/sblab/portel01/"\
    + "project/euni_5fc/intersect_histone_marks_cpgi/"

marks = ["H3K27ac", "H3K4me1", "H3K4me3"]
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


# Histone marks chipseq bigwig 
H3K27ac = 'H3K27acCHIP_ENCFF483XXE_ENCFF526YRA_mm9.bw'
H3K4me1 = 'H3K4me1CHIP_mm9.bw'
H3K4me3 = 'H3K4me3CHIP_mm9.bw'
hist_bw_files = OrderedDict((name, path_hist_bw + "/" + eval(name))
                              for name in ['H3K27ac', 'H3K4me1', 'H3K4me3'])

conditions = {"WT_histchipseq": bed_marks_WT, "KO_histchipseq": bed_marks_KO, "WT_KO_histchipseq": bed_marks_WT_KO}
subtype="histchipseq"


# Start mark iteration
for condition, dict_condition in conditions.items():
    for mark, mark_bed_file in dict_condition.items():

        # using templated printing options
        # mark_bed_file should allow for a list as well

        options_scale_reg = {"reference": " ".join(mark_bed_file),
                             "score_list": " ".join(hist_bw_files.values()),
                             "out_matrix": subtype + "_" + mark + "_" + condition +  "_scaleregions.mat.gz",
                             "start_label": "start",
                             "end_label": "end",
                             "region_len": 500,
                             "len_before": 1000,
                             "len_after": 1000
                             }

        options_plotprofile = {"matrix_file": subtype + "_" + mark + "_" + condition + "_scaleregions.mat.gz",
                               "out_plot": subtype + "_" + mark + "_" + condition +  "_scaleregions.pdf",
                               "out_dat": subtype + "_" + mark + "_" + condition + "_scaleregions.dat",
                               "plot_type": "se",
                               "colors": " ".join(["red", "blue", "green"]),
                               "start_label": "start",
                               "end_label": "end",
                               "samples_label": " ".join(hist_bw_files.keys()),
                               "regions_label": " ".join(['"' + mark + '"',
                                                          '"' + mark + '_no5fC"']),
                               "y_label": "RPKM"
                               }

        scale_reg = scaleregions(options_scale_reg)
        plot_reg = plotprofile(options_plotprofile)
        make_dir(condition)
        with cd(condition):
            make_dir(mark)
            with cd(mark):
                with open("submit_" + mark + ".sh", "w") as out:
                    cmd = scale_reg.cmd() + " && \ " + plot_reg.cmd()
                    out.write(cmd)
                    slurm.sbatch(cmd, J=mark, log=mark + ".log")
