# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from emtolib import walk_emtodirs
from emtolib.input_files import EmtoKgrnFile
from emtolib.config import read_config, update_slurm_settings, update_emto_paths

CONFIG = read_config()

SWS_NB = 3.071

PARAMS = {
    "niter": 200,
    "afm": "P",
    "ibz": 3,
    "nky": 25,
    "nz1": 32,
    "nz2": 8,
    "nz3": 8,
    "nres": 4,
    "nzd": 5000,
    "depth": 1.0,
    "imagz": 0.001,
    "eps": 0.2,
    "amix": 0.1,
    "efmix": 0.1,
    "tole": 1e-7,
    "tolef": 1e-7,
    "tolcpa": 1e-6,
    "sws": SWS_NB,
    "efgs": 0.1,
    "hx": 0.3,
    "nx": 7,
}


def diff_folders(root, ignore=("jobnam",)):
    folders = list(walk_emtodirs(root))
    folder1 = folders[0]
    dat1 = folder1.dat
    diffs = dict()
    for folder in folders[1:]:
        diff = dat1.param_diff(folder.dat)
        for k in ignore:
            if k in diff:
                diff.pop(k)
        if diff:
            diffs[folder.root] = diff

    for k, diff in diffs.items():
        print(k, diff)


def main():
    root = Path("app") / "Nb-Ta" / "CPAf"
    executable = str(CONFIG["emto"]["executable2"])

    # diff_folders(root, ignore=("jobnam", ))

    for folder in walk_emtodirs(root):
        print("Updating", folder.root)
        dat: EmtoKgrnFile = folder.dat
        slurm = folder.get_slurm()

        # Update KGRN input file
        update_emto_paths(dat, CONFIG, kstr="bcc.tfh", bmdl="bcc.mdl")
        dat.update(PARAMS)
        dat.strt = "A"
        dat.nky = 75
        c_nb = dat.get_concentration("Nb")

        # Update SLURM script
        update_slurm_settings(slurm, CONFIG, executable, dat.jobnam)
        slurm.jobname = f"NbTa-{int(100 * c_nb):02d}"

        # dat.dump()
        # slurm.dump()


if __name__ == "__main__":
    main()
