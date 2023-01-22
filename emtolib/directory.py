# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import os
import re

RE_COMP = re.compile(r"(\w+?)(\d+)")


class EmtoDirectory:
    def __init__(self, path):
        self.root = path

    @property
    def dat_path(self):
        paths = list()
        for name in os.listdir(self.root):
            if name.endswith(".dat"):
                paths.append(os.path.join(self.root, name))
        assert len(paths) == 1
        return paths[0]

    @property
    def prn_path(self):
        paths = list()
        for name in os.listdir(self.root):
            if name.endswith(".prn"):
                paths.append(os.path.join(self.root, name))
        if len(paths) > 1:
            paths.sort(key=lambda x: len(x))
        # assert len(paths) == 1
        return paths[0]

    @property
    def dos_paths(self):
        paths = list()
        for name in os.listdir(self.root):
            if name.endswith(".dos"):
                paths.append(os.path.join(self.root, name))
        paths.sort(key=lambda x: len(x))
        return paths

    # def dat_file(self):
    #     return EmtoDatFile(self.dat_path)

    # def prn_file(self):
    #     return EmtoPrnFile(self.prn_path)

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.root})>"


def parse_dirname(dirname, atoms):
    dirname = os.path.split(dirname)[1]
    atoms_lower = [at.lower() for at in atoms]
    match = RE_COMP.findall(dirname)
    concs = [0] * len(atoms)
    if match:
        for at, c in match:
            i = atoms_lower.index(at.lower())
            concs[i] = int(c)
    else:
        idx = atoms_lower.index(dirname.lower())
        concs[idx] = 100
    return tuple(concs)


def get_emto_dirnames(root, atoms):
    paths = dict()
    for name in os.listdir(root):
        path = os.path.join(root, name)
        try:
            cc = parse_dirname(name, atoms)
            paths[cc] = path
        except ValueError:
            pass
    return paths
