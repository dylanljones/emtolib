# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import re
import shutil
import warnings
from pathlib import Path
from typing import Union
from .files import KgrnFile, BmdlFile, PrnFile, DosFile, SlurmScript

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

        self._dat: Union[KgrnFile, None] = None
        self._prn: Union[PrnFile, None] = None
        self._dos: Union[DosFile, None] = None
        self._slurm: Union[SlurmScript, None] = None

    @property
    def dat(self):
        if self._dat is None:
            try:
                return self.get_dat()
            except FileNotFoundError:
                return None
        return self._dat

    @property
    def prn(self):
        if self._prn is None:
            try:
                return self.get_prn()
            except FileNotFoundError:
                return None
        return self._dat

    @property
    def dos(self):
        if self._dos is None:
            try:
                return self.get_dos()
            except FileNotFoundError:
                return None
        return self._dos

    @property
    def slurm(self):
        if self._slurm is None:
            try:
                return self.get_slurm()
            except FileNotFoundError:
                return None
        return self._slurm

    def get_dat_path(self):
        return find_input_file(self.path)

    def get_dat(self, path=""):
        if not path:
            path = self.get_dat_path()
        dat = KgrnFile(path)
        self._dat = dat
        return dat

    def get_input_path(self):
        warnings.warn(
            "get_input_path() is deprecated, use get_dat_path() instead",
            DeprecationWarning,
        )
        return self.get_dat_path()

    def get_input(self, path=""):
        warnings.warn(
            "get_input() is deprecated, use get_dat() instead", DeprecationWarning
        )
        return self.get_dat(path)

    def get_dos_path(self, name=""):
        if not name:
            name = self.dat.jobnam
        return self.path / f"{name}.dos"

    def get_dos(self, name=""):
        path = self.get_dos_path(name)
        dos = DosFile(path)
        self._dos = dos
        return dos

    def get_prn_path(self, name=""):
        if not name:
            name = self.dat.jobnam
        return self.path / f"{name}.prn"

    def get_prn(self, name=""):
        path = self.get_prn_path(name)
        prn = PrnFile(path)
        self._prn = prn
        return prn

    def get_slurm_path(self, name=""):
        name = name or "run_emto"
        return self.path / name

    def get_slurm(self, name=""):
        path = self.get_slurm_path(name)
        slurm = SlurmScript(path)
        self._slurm = slurm
        return slurm

    def get_slurm_out_paths(self):
        paths = list()
        for path in self.path.iterdir():
            if path.is_file():
                if path.name.startswith("slurm") and path.suffix == ".out":
                    paths.append(path)
        return paths

    def get_mdl_path(self, base=""):
        path = Path(self.dat.for004).with_suffix(".dat")
        if base:
            path = Path(base) / path.name
        return path

    def get_mdl(self, base=""):
        path = self.get_mdl_path(base)
        return BmdlFile(path)

    def move(self, dst):
        dst = Path(dst)
        if dst.exists():
            raise FileExistsError(f"Destination {dst} already exists!")
        shutil.move(self.path, dst)
        self.path = dst

    def rename(self, name: str):
        self.move(self.path.parent / name)

    def copy(self, dst, exist_ok=False):
        dst = Path(dst)
        if dst.exists():
            if not exist_ok:
                raise FileExistsError(f"Destination {dst} already exists!")
        else:
            shutil.copytree(self.path, dst)
        folder = self.__class__(dst)
        return folder

    def mkdirs(self):
        for name in self.dat.aux_dirs():
            path = self.path / name
            path.mkdir(parents=True, exist_ok=True)

    def clear(self, slurm=True, prn=True, dos=True, aux=True):
        dat = self.dat
        if slurm:
            for path in self.get_slurm_out_paths():
                path.unlink(missing_ok=True)
        if prn:
            path = self.get_prn_path()
            path.unlink(missing_ok=True)
        if dos:
            path = self.get_dos_path()
            path.unlink(missing_ok=True)
        if aux:
            for name in dat.aux_dirs():
                if name.startswith("/"):
                    name = name[1:]
                path = self.path / name
                if path.is_dir() and path.exists():
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
