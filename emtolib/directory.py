# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import os
import re
import shutil
from .input_files import EmtoKgrnFile
from .ouput_files import EmtoPrnFile, EmtoDosFile

RE_COMP = re.compile(r"(\w+?)(\d+)")


def searchdir(root_dir, pattern, recursive=False):
    pattern = pattern.replace(".", r"\.").replace("*", ".+")
    regex = re.compile(pattern)
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            path = os.path.join(root, file)
            relpath = os.path.relpath(path, root_dir)
            match = regex.match(relpath)
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

    def join(self, *relpath):
        return os.path.join(self.root, *relpath)

    def relpath(self, path):
        return os.path.relpath(path, self.root)

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
        self.dat = None
        try:
            self.dat = self.get_input()
        except FileNotFoundError:
            pass

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

    def get_dos_path(self, name=""):
        if not name:
            name = self.dat.jobname
        path = os.path.join(self.root, name + ".dos")
        return path

    def get_dos(self, name=""):
        path = self.get_dos_path(name)
        return EmtoDosFile(path)

    def get_prn_path(self, name=""):
        if not name:
            name = self.dat.jobname
        path = os.path.join(self.root, name + ".prn")
        return path

    def get_prn(self, name=""):
        path = self.get_prn_path(name)
        return EmtoPrnFile(path)


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
