# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-06-21

import shutil
import logging
from typing import Union
from pathlib import Path

import numpy as np
from scipy import constants as cc

from .files import (
    DosFile,
    PrnFile,
    BmdlFile,
    DmftFile,
    KgrnFile,
    SlurmScript,
    read_sigma_z,
    read_sigma_iw,
    read_fermi_surface_file,
)
from .errors import KGRNReadError
from .files.dos_file import DosLmsFile

logger = logging.getLogger(__name__)

ry2ev = cc.value("Rydberg constant times hc in eV")  # 13.605693122994


def find_kgrn_file(folder: Union[Path, str]) -> Path:
    folder = Path(folder)
    # get *.dat files
    for file in folder.glob("*.dat"):
        if file.is_file() and not file.name.startswith("dmft"):
            # Check contents of file
            with open(file, "r") as fh:
                fkey = fh.read(4)
                if fkey in ("KGRN", "DMFT"):
                    return file
    raise FileNotFoundError(f"No input file found in {folder}!")


class EmtoDirectory:
    """Class to handle EMTO simulation directories."""

    def __init__(self, *path, missing_ok=True):
        self.path = Path(*path)
        if not missing_ok and not self.path.exists():
            raise FileNotFoundError(f"Directory {self.path} does not exist!")

        self._dat: Union[KgrnFile, None] = None
        self._prn: Union[PrnFile, None] = None
        self._dos: Union[DosFile, None] = None
        self._dos_lms: Union[DosLmsFile, None] = None
        self._slurm: Union[SlurmScript, None] = None
        self._dmft: Union[DmftFile, None] = None

    @property
    def name(self):
        return self.path.name

    @property
    def parent(self):
        return self.path.parent

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
        return self._prn

    @property
    def dos(self):
        if self._dos is None:
            try:
                return self.get_dos()
            except FileNotFoundError:
                return None
        return self._dos

    @property
    def dos_lms(self):
        if self._dos_lms is None:
            try:
                return self.get_dos_lms()
            except FileNotFoundError:
                return None
        return self._dos_lms

    @property
    def slurm(self):
        if self._slurm is None:
            try:
                return self.get_slurm()
            except FileNotFoundError:
                return None
        return self._slurm

    @property
    def dmft(self):
        if self._dmft is None:
            try:
                return self.get_dmft()
            except FileNotFoundError:
                return None
        return self._dmft

    def has_input_file(self, allow_ini=True):
        try:
            if self.dat is not None:
                return True
            elif allow_ini and (self.path / "input.ini").exists():
                return True
        except KGRNReadError as e:
            logger.exception(e)
        return False

    def converged(self):
        if self.prn is None:
            return False
        return self.prn.converged

    def get_dat_path(self):
        if not self.path.exists():
            raise FileNotFoundError(f"Directory {self.path} does not exist!")
        return find_kgrn_file(self.path)

    def get_jobnam(self):
        if not self.path.exists():
            raise FileNotFoundError(f"Directory {self.path} does not exist!")
        if self.dat is None:
            raise FileNotFoundError(f"No input file found in {self.path}!")
        return self.dat.jobnam

    def get_dos_path(self, name=""):
        if not name:
            name = self.get_jobnam()
        return self.path / f"{name}.dos"

    def get_dos_lms_path(self, name=""):
        if not name:
            name = self.get_jobnam()
        return self.path / f"{name}.lms"

    def get_fes_path(self, name=""):
        if not name:
            name = self.get_jobnam()
        return self.path / f"{name}.fes"

    def get_prn_path(self, name=""):
        if not name:
            name = self.get_jobnam()
        return self.path / f"{name}.prn"

    def get_mdl_path(self, base=""):
        path = Path(self.dat.for004).with_suffix(".dat")
        if base:
            path = Path(base) / path.name
        return path

    def get_slurm_path(self, name=""):
        name = name or "run_emto"
        return self.path / name

    def get_slurm_out_paths(self):
        paths = list()
        for path in self.path.iterdir():
            if path.is_file():
                if path.name.startswith("slurm") and path.suffix == ".out":
                    paths.append(path)
        paths.sort(key=lambda x: x.name.replace("slurm-", "").replace(".out", ""))
        return paths

    def get_dmft_path(self, name=""):
        if not name:
            name = "dmft.dat"
        path = self.path / name
        return path

    def get_dat(self, path="", reload=False):
        if not reload and self._dat is not None:
            return self._dat
        if not path:
            path = self.get_dat_path()
        else:
            path = self.path / path
        dat = KgrnFile(path)
        self._dat = dat
        return dat

    def get_dos(self, name="", reload=False):
        if not reload and self._dos is not None:
            return self._dos
        path = self.get_dos_path(name)
        dos = DosFile(path)
        self._dos = dos
        return dos

    def get_dos_lms(self, name="", reload=False):
        if not reload and self._dos_lms is not None:
            return self._dos_lms
        path = self.get_dos_lms_path(name)
        dos = DosLmsFile(path)
        self._dos_lms = dos
        return dos

    def get_prn(self, name="", reload=False):
        if not reload and self._prn is not None:
            return self._prn
        path = self.get_prn_path(name)
        prn = PrnFile(path)
        self._prn = prn
        return prn

    def get_slurm(self, name="", reload=False):
        if not reload and self._slurm is not None:
            return self._slurm
        path = self.get_slurm_path(name)
        slurm = SlurmScript(path)
        self._slurm = slurm
        return slurm

    def get_mdl(self, base=""):
        path = self.get_mdl_path(base)
        return BmdlFile(path)

    def get_dmft(self, name=""):
        path = self.get_dmft_path(name)
        return DmftFile(path)

    def exists(self):
        return self.path.exists()

    def mkdir(self, parents=True, exist_ok=True):
        self.path.mkdir(parents=parents, exist_ok=exist_ok)

    def rmdir(self, missing_ok=True):
        if not missing_ok and not self.path.exists():
            raise FileNotFoundError(f"Directory {self.path} does not exist!")
        shutil.rmtree(self.path)

    def iter_dir(self):
        return self.path.iterdir()

    def glob(self, pattern="*", recursive=False):
        if recursive:
            return self.path.rglob(pattern)
        else:
            return self.path.glob(pattern)

    def relpath(self, path):
        return self.path.relative_to(path)

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

    def mkdirs(self, keep=False):
        for name in self.dat.aux_dirs():
            path = self.path / name
            path.mkdir(parents=True, exist_ok=True)
            if keep:
                keep = path / ".keep"
                keep.touch(exist_ok=True)

    def clear(self, slurm=True, out=True, aux=True):
        dat = self.dat
        if slurm:
            for path in self.get_slurm_out_paths():
                path.unlink(missing_ok=True)
        if out:
            jobname = dat.jobnam
            # prn file
            file = self.get_prn_path()
            file.unlink(missing_ok=True)
            # dos files
            file = self.get_dos_path()
            file.unlink(missing_ok=True)
            file = self.get_dos_lms_path()
            file.unlink(missing_ok=True)
            # fermi surface
            file = self.get_fes_path()
            file.unlink(missing_ok=True)
            # other out files
            file = self.path / f"{jobname}.phi"
            file.unlink(missing_ok=True)
            file = self.path / f"{jobname}.epm"
            file.unlink(missing_ok=True)
            # fort files
            for file in self.path.iterdir():
                if file.name.startswith("fort"):
                    file.unlink(missing_ok=True)
                elif file.name == "Sig":
                    file.unlink(missing_ok=True)
        if aux:
            for name in dat.aux_dirs():
                if name.startswith("/"):
                    name = name[1:]
                path = self.path / name
                if path.is_dir() and path.exists():
                    for file in path.iterdir():
                        file.unlink(missing_ok=True)

    def clear_slurm(self):
        for path in self.get_slurm_out_paths():
            path.unlink(missing_ok=True)

    def __truediv__(self, other):
        new_path = self.path / other
        if new_path.is_dir():
            return self.__class__(new_path)
        else:
            return new_path

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.path})>"

    def __str__(self):
        return str(self.path)

    def convert_dmft_to_kgrn(self, dst=None):
        assert self.dat.is_dmft
        folder = self if dst is None else self.copy(dst, exist_ok=False)
        dat = folder.dat

        # Create dmft.dat file from KGRN-DMFT parameters
        dmft = DmftFile(folder / "dmft.dat")
        dmft.for007 = dat.for007
        dmft.nom = dat.nom
        dmft.ttt = dat.ttt
        dmft.solver = dat.solver
        dmft.dc = dat.dc
        dmft.smix = dat.smix
        dmft.nomi = dat.nomi

        # Save files
        dat.force_dmft(False)
        dat.dump()
        dmft.dump()

        return folder

    def convert_kgrn_to_dmft(self, dst=None):
        assert not self.dat.is_dmft
        folder = self if dst is None else self.copy(dst, exist_ok=False)
        dat = folder.dat
        dmft = folder.dmft

        # Update parameters of KGRN-DMFT from dmft.dat file
        dat.update(dmft.to_dict())

        # remove dmft.dat file
        dmft.path.unlink()
        folder._dmft = None

        # Save *.dat file
        dat.force_dmft(True)
        dat.dump()

    def get_aimag(self):
        num = 300
        path = self.path / f"fort.{num}"
        if not path.exists():
            raise FileNotFoundError(f"File {path.name} does not exist!")
        return np.loadtxt(path)

    def get_sigma_z(self, name="", unit="ry"):
        num = 99
        if not name:
            name = f"fort.{num}"
        path = self.path / name
        return read_sigma_z(path, unit)

    def get_sigma_iw(self, name="", unit="ry"):
        num = 400
        if not name:
            name = f"fort.{num}"
        path = self.path / name
        return read_sigma_iw(path, unit)

    def get_fermi_surface(self, name=""):
        path = self.get_fes_path(name)
        return read_fermi_surface_file(path)


def is_emtodir(path):
    """Checks if path is an EMTO directory containing at least the input *.dat file."""
    path = Path(path)
    if not path.is_dir():
        return False
    try:
        find_kgrn_file(path)
        return True
    except FileNotFoundError:
        return False


def walk_emtodirs(*paths, recursive=False, missing_dat_ok=False):
    """Walks through all emto directories in the given paths.

    Parameters
    ----------
    paths : str or Path
        The paths to walk through. Note that the paths can contain an EMTO
        directory itself.
    recursive : bool, optional
        If True, the paths will be walked recursively.
    missing_dat_ok : bool, optional
        If True, folders without a dat file will be included. This is equivalent to
        iterating over all directories.
    """
    paths = paths or (".",)  # Use current directory if no path is given
    for root in paths:
        root = Path(root)
        # Check if the given root path is an EMTO directory
        if root.is_dir():
            folder = EmtoDirectory(root)
            if folder.has_input_file() or missing_dat_ok:
                yield folder

        # Iterate (recursively) over all folders in the given root path
        iterator = root.rglob("*") if recursive else root.glob("*")
        for folder in iterator:
            if folder.is_dir():
                folder = EmtoDirectory(folder)
                if folder.has_input_file() or missing_dat_ok:
                    yield folder


def walk_ini_emtodirs(*paths, recursive=False):
    """Walks through all emto directories in the given paths.

    Parameters
    ----------
    paths : str or Path
        The paths to walk through. Note that the paths can contain an EMTO
        directory itself.
    recursive : bool, optional
        If True, the paths will be walked recursively.
    """
    paths = paths or (".",)  # Use current directory if no path is given
    for root in paths:
        root = Path(root)
        # Check if the given root path is an EMTO directory
        if root.is_dir():
            folder = EmtoDirectory(root)
            if len(list(folder.glob("*.ini"))):
                yield folder

        # Iterate (recursively) over all folders in the given root path
        iterator = root.rglob("*") if recursive else root.glob("*")
        for folder in iterator:
            if folder.is_dir():
                folder = EmtoDirectory(folder)
                if len(list(folder.glob("*.ini"))):
                    yield folder


def diff_emtodirs(*paths, recursive=False, exclude=None):
    """Find difference of the input files of all EMTO directories in root.

    Does not support diffing the atom configs, only the parameters!
    """
    # Find all parameters that change
    changing_params = set()
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    dat_base = folders[0].dat
    for folder in folders[1:]:
        diff = dat_base.param_diff(folder.dat, exclude)
        changing_params.update(diff.keys())

    # Exract values of changing parameters
    if not changing_params:
        return dict()
    diffs = dict()
    for folder in folders:
        diff = {k: folder.dat[k] for k in changing_params}
        diffs[folder.path] = diff

    return diffs
