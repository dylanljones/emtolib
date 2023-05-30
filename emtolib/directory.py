# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import re
import shutil
from pathlib import Path
from .input_files import EmtoKgrnFile
from .input_files.kgrn_dmft import EmtoKgrnFile as EmtoKgrnFileDMFT
from .ouput_files import EmtoPrnFile, EmtoDosFile
from typing import Union

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

    def __init__(self, path, dmft=False):
        self.root = Path(path)
        self.dat = None
        try:
            self.dat = self.get_input(dmft=dmft)
        except FileNotFoundError:
            pass

    def move(self, dst):
        dst = Path(dst)
        shutil.move(self.root, dst)
        self.root = dst

    def copy(self, dst):
        dst = Path(dst)
        shutil.copytree(self.root, dst)
        folder = self.__class__(dst)
        return folder

    def get_input_path(self):
        return find_input_file(self.root)

    def get_input(self, path="", dmft=False):
        if not path:
            path = self.get_input_path()
        if dmft:
            file = EmtoKgrnFileDMFT(path)
        else:
            file = EmtoKgrnFile(path)
        return file

    def get_dos_path(self, name=""):
        if not name:
            name = self.dat.jobname
        return self.root / f"{name}.dos"

    def get_dos(self, name=""):
        path = self.get_dos_path(name)
        return EmtoDosFile(path)

    def get_prn_path(self, name=""):
        if not name:
            name = self.dat.jobname
        return self.root / f"{name}.prn"

    def get_prn(self, name=""):
        path = self.get_prn_path(name)
        return EmtoPrnFile(path)

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.root})>"
