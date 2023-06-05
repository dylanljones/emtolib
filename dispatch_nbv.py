# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from emtolib import EmtoDirectory, walk_emtodirs
from emtolib.config import read_config

CONFIG = read_config()


def main():
    executable = str(CONFIG["emto"]["executable2"])

    root = Path("app") / "Nb-V" / "CPA2"
    for folder in walk_emtodirs(root):
        print("Updating", folder.root)
        dat = folder.dat
        dat.strt = "A"
        dat.nky = 57
        dat.dump()


if __name__ == "__main__":
    main()
