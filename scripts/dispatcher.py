# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-20

from emtolib import EmtoDirectory, Atom, elements, CONFIG
from emtolib.config import update_slurm_settings, update_emto_paths

ROOT = CONFIG["data_dir"]
EXEC_DMFT = CONFIG["emto"]["executable"]

PARAMS_NB = {
    "niter": 999,
    "afm": "P",
    "ibz": 3,
    "nky": 25,
    "nz1": 32,
    "nz2": 8,
    "nz3": 8,
    "nres": 4,
    "nzd": 2000,
    "depth": 1.0,
    "imagz": 0.001,
    "eps": 0.2,
    "amix": 0.1,
    "efmix": 0.1,
    "tole": 1e-7,
    "tolef": 1e-7,
    "tolcpa": 1e-6,
    "efgs": 0.1,
    "hx": 0.3,
    "nx": 7,
    "ncpa": 30,
}

PARAMS_V = PARAMS_NB.copy()


def main():
    at = "Nb"
    nl = 3
    root = ROOT / at

    folder = EmtoDirectory(root / f"NL{nl}")
    folder.mkdir(exist_ok=True)
    if folder.dat is None:
        dat = folder.get_dat(f"{at.lower()}.dat")
    else:
        dat = folder.dat
    slurm = folder.get_slurm()
    dmft = folder.get_dmft()
    print(dmft)

    dat.jobnam = at.lower()
    dat.atoms.clear()
    dat.add_atom(Atom(at))
    dat.update_mnta()

    update_emto_paths(dat, CONFIG, f"bcc_nl{nl}.tfh", bmdl="bcc.mdl")
    dat.update(PARAMS_NB)
    dat.strt = "A"
    dat.sws = elements[at]["sws"]

    update_slurm_settings(slurm, CONFIG)
    slurm.set_body(EXEC_DMFT, dat.path.name)
    slurm.jobname = f"{at}_NL{nl}"

    dat.dump()
    slurm.dump()
    folder.clear(keep=True)


if __name__ == "__main__":
    main()
