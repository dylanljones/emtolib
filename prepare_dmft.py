# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from emtolib.input_files.kgrn_dmft import EmtoKgrnFile
from emtolib.directory import find_input_file


def get_kgrn_files(root):
    root = Path(root)
    for folder in root.iterdir():
        if folder.is_dir():
            inp = find_input_file(folder)
            dat = EmtoKgrnFile(inp)
            yield dat


def update_kgrn(dat, ud):
    dat.strt = "A"
    dat.expan = "M"
    dat.for001_2 = "../../../EMTO5.8/kstr/smx/bcc-10.tfh"
    dat.niter = 200
    dat.nky = 7
    dat.ibz = 3
    dat.ibz2 = 1
    # dat.shf = 0
    dat.nz1 = 32
    dat.nz2 = 8
    dat.nz3 = 8
    dat.nres = 4
    dat.nzd = 3000
    dat.depth = 1.4
    dat.imagz = 0.01
    # dat.eps = 0.2
    dat.amix = 0.1
    dat.efmix = 1.0
    dat.efgs = 0.0
    dat.hx = 0.3
    dat.nx = 7
    dat.nz0 = 16

    for at in dat.atoms:
        if at.symbol == "Nb":
            at.u = [0, 0, 0, 0]
            at.j = [0, 0, 0, 0]
            at.valen[-4:] = 1
        elif at.symbol == "V":
            at.u = [0, 0, ud, 0]
            at.j = [0, 0, 0, 0]


def main():
    root = Path("")
    for dat in get_kgrn_files(root):
        print(dat.path)
        update_kgrn(dat, 2.0)
        dat.dump()


if __name__ == "__main__":
    main()
