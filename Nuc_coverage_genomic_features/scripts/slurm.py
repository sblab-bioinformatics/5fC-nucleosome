#! /usr/bin/env python3

import subprocess
import sys

def sbatch(cmd, J="fortytwo", log="fortywo.log", mem=65536 ):
    if len(cmd) == 0:
        print("Error, command can not be empy")
        sys.error
    sbatch_cmd = 'sbatch -J ' + J + ' -o ' + log + ' --mem ' + str(mem) +\
            ' --wrap="' + cmd + '"'
    subprocess.call(sbatch_cmd, shell=True)
    return

