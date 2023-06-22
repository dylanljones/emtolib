# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import re
import shutil
from pathlib import Path
from typing import Union
from .input_files import EmtoKgrnFile
from .output_files import EmtoPrnFile, EmtoDosFile
from .slurm import SlurmScript

RE_COMP = re.compile(r"(\w+?)(\d+)")


def find_input_file(folder: Union[Path, str]) -> Path:
    folder = Path(folder)
    # get *.dat files
    for file in folder.glob("*.dat"):
        if file.is_file() and not folder.name.startswith("dmft"):
            # Check contents of file
            text = file.read_text().strip()
            if text.startswith("KGRN"):
                return file
    raise FileNotFoundError(f"No input file found in {folder}!")


class EmtoDirectory:
    """Class to handle EMTO simulation directories."""

    def __init__(self, path):
        self.path = Path(path)
        self.dat: Union[EmtoKgrnFile, None] = None
        try:
            self.dat = self.get_input()
        except FileNotFoundError:
            self.dat = None

    def move(self, dst):
        dst = Path(dst)
        shutil.move(self.path, dst)
        self.path = dst

    def copy(self, dst):
        dst = Path(dst)
        shutil.copytree(self.path, dst)
        folder = self.__class__(dst)
        return folder

    def get_input_path(self):
        return find_input_file(self.path)

    def get_input(self, path=""):
        if not path:
            path = self.get_input_path()
        file = EmtoKgrnFile(path)
        return file

    def get_dos_path(self, name=""):
        if not name:
            name = self.dat.jobnam
        return self.path / f"{name}.dos"

    def get_dos(self, name=""):
        path = self.get_dos_path(name)
        return EmtoDosFile(path)

    def get_prn_path(self, name=""):
        if not name:
            name = self.dat.jobnam
        return self.path / f"{name}.prn"

    def get_slurm_out_paths(self):
        paths = list()
        for path in self.path.iterdir():
            if path.is_file():
                if path.name.startswith("slurm") and path.suffix == ".out":
                    paths.append(path)
        return paths

    def get_prn(self, name=""):
        path = self.get_prn_path(name)
        return EmtoPrnFile(path)

    def get_slurm(self, name="run_emto"):
        path = self.path / name
        return SlurmScript(path)

    def mkdirs(self):
        for name in self.dat.aux_dirs():
            path = self.path / name
            path.mkdir(parents=True, exist_ok=True)

    def clear(self, slurm=True, prn=True, dos=True, aux=True):
        dat = self.dat
        if slurm:
            for path in self.get_slurm_out_paths():
                path.unlink()
        if prn:
            path = self.get_prn_path()
            path.unlink()
        if dos:
            path = self.get_dos_path()
            path.unlink()
        if aux:
            for name in dat.aux_dirs():
                path = self.path / name
                if path.exists():
                    shutil.rmtree(path)
            self.mkdirs()

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.path})>"


def walk_emtodirs(root):
    root = Path(root)
    for folder in root.glob("*"):
        if not folder.is_dir():
            continue
        folder = EmtoDirectory(folder)
        if folder.dat is None:
            continue
        yield folder
