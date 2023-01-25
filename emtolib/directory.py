# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import os
import re
import shutil
from .input_files import EmtoKgrnFile

RE_COMP = re.compile(r"(\w+?)(\d+)")


def searchdir(root_dir, pattern, recursive=False):
    pattern = pattern.replace(".", r"\.").replace("*", ".+")
    regex = re.compile(pattern)
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            path = os.path.join(root, file)
            match = regex.match(path)
            if match:
                yield path
        if not recursive:
            break


class Directory:

    def __init__(self, path):
        self.root = path

    def move(self, dst):
        shutil.move(self.root, dst)
        self.root = dst

    def copy(self, dst):
        shutil.copytree(self.root, dst)
        folder = self.__class__(dst)
        return folder

    def listdir(self):
        return os.listdir(self.root)

    def searchdir(self, pattern, recursive=False):
        return list(searchdir(self.root, pattern, recursive))

    def listfiles(self):
        for name in os.listdir(self.root):
            path = os.path.join(self.root, name)
            if os.path.isfile(path):
                yield path

    def listdirs(self):
        for name in os.listdir(self.root):
            path = os.path.join(self.root, name)
            if os.path.isdir(path):
                yield path

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.root})>"


# =========================================================================


class EmtoDirectory(Directory):

    def __init__(self, path):
        super().__init__(path)

    def get_input_path(self):
        dat_paths = list(self.searchdir("*.dat"))
        for p in dat_paths[:]:
            if os.path.split(p)[1].startswith("dmft"):
                dat_paths.remove(p)
        if len(dat_paths) != 1:
            raise ValueError("Could not identify input '*.dat' file")
        return dat_paths[0]

    def get_input(self, path=""):
        if not path:
            path = self.get_input_path()
        return EmtoKgrnFile(path)


class EmtoDirectoryOld:
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
