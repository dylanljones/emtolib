# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from emtolib import walk_emtodirs, EmtoDirectory
from emtolib.input_files import EmtoKgrnFile
from emtolib.config import read_config
from emtolib.config import update_slurm_settings, update_emto_paths

CONFIG = read_config()


def update(dat: EmtoKgrnFile, slurm, executable):
    update_emto_paths(dat, CONFIG, kstr="bcc.tfh", bmdl="bcc.mdl")
    dat.strt = "A"
    dat.nzd = 5000
    dat.nky = 17
    dat.depth = 1.0
    dat.imagz = 0.001
    dat.eps = 0.2
    dat.amix = 0.01
    dat.efmix = 1.0
    dat.efgs = 0.1
    update_slurm_settings(slurm, CONFIG, executable, dat.jobname)


def main():
    root = Path("app") / "Nb-Ta" / "CPA2"
    executable = str(CONFIG["emto"]["executable2"])

    for folder in walk_emtodirs(root):
        dat: EmtoKgrnFile = folder.dat
        dat.strt = "N"
        dat.nky = 75
        dat.hx = 0.3
        dat.nx = 7
        dat.dump()

        # slurm = folder.get_slurm()
        # slurm.dump()


if __name__ == "__main__":
    main()
