# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-04

import pprint
from pathlib import Path
from emtolib import walk_emtodirs, EmtoKgrnFile, elements, CONFIG
from emtolib.config import update_slurm_settings, update_emto_paths
from emtolib.files.make_file import generate_makefile

DEFAULT_PARAMS = {
    "jobnam": "nb",
    "strt": "A",
    "msgl": 1,
    "expan": "S",
    "fcd": "Y",
    "func": "SCA",
    "for001": "~/EMTO/EMTO5.8/kstr/smx/bcc.tfh",
    "for001_2": "",
    "for004": "~/EMTO/EMTO5.8/bmdl/mdl/bcc.mdl",
    "dir002": "pot/",
    "dir003": "pot/",
    "dir006": "",
    "dir009": "pot/",
    "dir010": "chd/",
    "dir011": "/tmp/",
    "niter": 500,
    "nlin": 31,
    "nprn": 0,
    "ncpa": 30,
    "nt": 1,
    # "mnta": 2,
    "mode": "3D",
    "frc": "N",
    "dos": "D",
    "ops": "N",
    "afm": "P",
    "crt": "M",
    "lmaxh": 8,
    "lmaxt": 4,
    "nfi": 31,
    "fixg": 2,
    "shf": 0,
    "sofc": "N",
    "kmsh": "G",
    "ibz": 3,
    "nkx": 0,
    "nky": 75,
    "nkz": 0,
    "fbz": "N",
    "kmsh2": "G",
    "ibz2": 1,
    "nkx2": 4,
    "nky2": 0,
    "nkz2": 51,
    "zmsh": "C",
    "nz1": 32,
    "nz2": 8,
    "nz3": 8,
    "nres": 4,
    "nzd": 5000,
    "depth": 1.0,
    "imagz": 0.01,
    "eps": 0.2,
    "elim": -1.0,
    "amix": 0.1,
    "efmix": 1.0,
    "vmtz": 0.0,
    "mmom": 0.0,
    "tole": 1e-7,
    "tolef": 1e-7,
    "tolcpa": 1e-6,
    "tfermi": 500.0,
    "sws": 2.8388,
    "nsws": 1,
    "dsws": 0.05,
    "alpcpa": 0.602,
    # ----------------
    "efgs": 0.0,
    "hx": 0.3,
    "nx": 9,
    "nz0": 16,
    "stmp": "N",
    # ----------------
    "iex": 4,
    "np": 251,
    "nes": 15,
    "dirac_niter": 100,
    "iwat": 0,
    "nprna": 0,
    "vmix": 0.3,
    "rwat": 3.5,
    "rmax": 20.0,
    "dx": 0.03,
    "dr1": 0.002,
    "test": 1e-12,
    "teste": 1e-12,
    "testy": 1e-12,
    "testv": 1e-12,
}


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
    "efgs": 0.1,
    "hx": 0.3,
    "nx": 7,
    "ncpa": 30,
}


def mix_sws(c, sws_c, sws_1mc):
    return c * sws_c + (1 - c) * sws_1mc


def param_diff(dat, data2):
    data1 = dat.to_dict()
    for k, v in data1.items():
        if k in data2 and v != data2[k]:
            print(f"{str(dat.path):<30} {k:<10} {v:<10} {data2[k]}")


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
    at1, at2 = "Nb", "V"
    nl = 4
    executable = str(CONFIG["emto"]["executable"])
    root = Path("app") / f"{at1}-{at2}" / f"CPA_{nl}_10"

    # base = EmtoDirectory(root / "Nb20")
    # data = base.dat.to_dict()
    # data = {k: data[k] for k in base.dat.order}
    # pprint.pprint(data, indent=4, sort_dicts=False)

    for folder in walk_emtodirs(root):
        print("Updating", folder.path)
        folder.clear(aux=False)
        dat: EmtoKgrnFile = folder.dat
        slurm = folder.get_slurm()
        mnta = len(dat.atoms)

        # Update KGRN input file
        update_emto_paths(dat, CONFIG, kstr=f"bcc_{nl}_10.tfh", bmdl="bcc.mdl")
        dat.update(PARAMS)
        dat.strt = "A"
        dat.nky = 57
        dat.mnta = mnta
        dat.for001_2 = ""

        c_nb = dat.get_concentration(at1)
        dat.sws = mix_sws(c_nb, elements.Nb.sws, elements.V.sws)

        # Update SLURM script
        update_slurm_settings(slurm, CONFIG, executable, dat.jobnam)
        slurm.jobname = f"{at1}{at2}-{int(100 * c_nb):02d}_{nl}_10"
        slurm.mem = "1gb"

        # --------------------------------------------------------------
        # param_diff(dat, DEFAULT_PARAMS)
        folder.mkdirs()
        dat.set_header()
        dat.dump()
        slurm.dump()

    makefile = generate_makefile(root)
    makefile.dump()


if __name__ == "__main__":
    main()
