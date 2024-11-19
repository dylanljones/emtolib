# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-17

import io
import re
from typing import TextIO

import numpy as np
import pandas as pd
from scipy import constants

from ..common import EmtoFile
from ..errors import DOSReadError

RY2EV = constants.value("Rydberg constant times hc in eV")  # 13.605693122994

RE_DOS_TOTAL = re.compile(r"Total DOS and NOS and partial \(IT\) (?P<spin>.*)$")
RE_DOS_SUBLAT = re.compile("Sublattice (?P<idx>.*) Atom (?P<atom>.*) spin (?P<spin>.*)")

UNITS = "ry", "ev"


def _search_sections_tnos(fp):
    fp.seek(0, io.SEEK_END)
    size = fp.tell()
    fp.seek(0)

    sec, typ = "", ""
    keys, sections = list(), list()
    while fp.tell() < size - 1:
        pos = fp.tell()
        line = fp.readline().strip(" \t\n#")
        end_pos = fp.tell()

        if line.startswith("Total"):
            sec, typ = line, "TDOS"
            keys.append((sec, typ))
            sections.append([end_pos, -1])
            if len(sections) > 1:
                sections[-2][1] = pos

        elif line.startswith("Sublatt"):
            sec, typ = line, "PDOS"
            keys.append((sec, typ))
            sections.append([end_pos, -1])
            if len(sections) > 1:
                sections[-2][1] = pos

        elif line.startswith("TNOS"):
            typ = "TNOS"
            keys.append((sec, typ))
            sections.append([end_pos, -1])
            if len(sections) > 1:
                sections[-2][1] = pos
    sections[-1][1] = size
    return dict(zip(keys, sections))


def _search_sections(fp: TextIO):
    fp.seek(0, io.SEEK_END)
    size = fp.tell()
    fp.seek(0)

    keys, sections = list(), list()
    while fp.tell() < size - 1:
        pos = fp.tell()
        line = fp.readline().strip(" \t\n#")
        end_pos = fp.tell()

        if line.startswith("Total"):
            sec = line
            keys.append(sec)
            sections.append([end_pos, -1])
            if len(sections) > 1:
                sections[-2][1] = pos

        elif line.startswith("Sublatt"):
            sec = line
            keys.append(sec)
            sections.append([end_pos, -1])
            if len(sections) > 1:
                sections[-2][1] = pos
    sections[-1][1] = size
    return dict(zip(keys, sections))


def _parse_columns(lines: list) -> str:
    """Parses the column line at the top of a section.

    Removes column line if it exists, so lines start at actual data
    """
    # check if columns exist
    header_exists = False
    for line in lines:
        if line:
            header_exists = any(c.isalpha() for c in line)
            break
    # Parse columns line
    columns = ""
    if header_exists:
        columns = lines.pop(0)
        while not columns:
            columns = lines.pop(0)
    return columns


def _parse_data(lines: list) -> tuple:
    dos_data, tnos_data = list(), list()
    contains_tnos = False
    i = 0
    for i, line in enumerate(lines):
        if line:
            if line.startswith("TNOS"):
                contains_tnos = True
                break
            values = list()
            for value in line.split():
                try:
                    values.append(float(value))
                except ValueError:
                    values.append(np.nan)
            dos_data.append(values)
    if contains_tnos:
        for line in lines[i + 1 :]:
            if line:
                values = [float(x) for x in line.split()]
                tnos_data.append(values)
    return dos_data, tnos_data


def read_dos(fp: TextIO):
    """Reads the contents of a EMTO DOS file.

    Parameters
    ----------
    fp : TextIO
        The opened file pointer.
    """
    sections = _search_sections(fp)
    tdos_data = list()
    pdos_data = list()
    tnos_data = list()
    for key, section in sections.items():
        start, end = section
        # Read data of section
        fp.seek(start)
        lines = [line.strip(" \t#") for line in fp.read(end - start).splitlines()]
        if key.startswith("Tot"):
            # Total DOS section
            match = RE_DOS_TOTAL.search(key)
            spin = match.group("spin")
            if spin.startswith("DOS"):
                spin = spin[3:]
            # spin = match.group("spin").removeprefix("DOS")

            _parse_columns(lines)
            columns = ["E", "Total", "TNOS", "Partial"]
            tdos, _ = _parse_data(lines)
            for values in tdos:
                row = {k: v for k, v in zip(columns, values)}
                row.update({"Spin": spin})
                tdos_data.append(row)

        elif key.startswith("Sub"):
            # Sublattice DOS section
            match = RE_DOS_SUBLAT.search(key)
            sublatt = int(match.group("idx"))
            atom = match.group("atom").strip()
            spin = match.group("spin").strip()
            index = {"Sublatt": sublatt, "Atom": atom, "Spin": spin}

            column_str = _parse_columns(lines)
            columns = column_str.split()
            pdos, tnos = _parse_data(lines)
            for values in pdos:
                row = {k: v for k, v in zip(columns, values)}
                row.update(index)
                pdos_data.append(row)
            for values in tnos:
                row = {k: v for k, v in zip(columns, values)}
                row.update(index)
                tnos_data.append(row)

    tdos = pd.DataFrame(tdos_data)
    tdos.set_index(["Spin"], inplace=True)
    pdos = pd.DataFrame(pdos_data)
    pdos.set_index(["Sublatt", "Atom", "Spin"], inplace=True)
    tnos = pd.DataFrame(tnos_data)
    tnos.set_index(["Sublatt", "Atom", "Spin"], inplace=True)

    return tdos, pdos, tnos


def load_dos(file_path: str):
    """Open and read the contents of a EMTO DOS file.

    Parameters
    ----------
    file_path : str
        The path of the EMTO DOS file.
    """
    with open(file_path, "r") as fp:
        return read_dos(fp)


def _index_dataframe(df, args, drop):
    if drop:
        i = 0
        for arg in args:
            if arg is not None:
                try:
                    df = df.xs(arg, level=i, drop_level=True)
                except TypeError:
                    df = df.loc[arg]
                    df.reset_index(drop=True, inplace=True)
            else:
                i += 1
    else:
        for i, arg in enumerate(args):
            if arg is not None:
                df = df.xs(arg, level=i, drop_level=False)
    return df


def translate_spin(spin):
    if spin is not None:
        spin = spin.upper()
        if spin == "DN":
            spin = "DOWN"
    return spin


def dos_ry2ev(energy, dos):
    return energy * RY2EV, dos / RY2EV



class DosFile(EmtoFile):
    extension = ".dos"

    def __init__(self, *path):
        super().__init__(*path)
        try:
            self.tdos, self.pdos, self.tnos = load_dos(self.path)
        except Exception as e:
            raise DOSReadError(f"Could not read DOS file: {self.path}") from e

    def get_tdos(self, spin=None, drop=True):
        spin = translate_spin(spin)
        return _index_dataframe(self.tdos, [spin], drop)

    def get_pdos(self, sublatt=None, atom=None, spin=None, drop=True):
        spin = translate_spin(spin)
        return _index_dataframe(self.pdos, [sublatt, atom, spin], drop)

    def get_tnos(self, sublatt=None, atom=None, spin=None, drop=True):
        spin = translate_spin(spin)
        return _index_dataframe(self.tnos, [sublatt, atom, spin], drop)

    def get_total_dos(self, spin=None, unit="ry"):
        unit = unit.lower()
        if unit not in UNITS:
            raise ValueError(f"Invalid unit: {unit}")

        df = self.get_tdos(spin, drop=True)
        series = df.groupby("E")["Total"].sum()
        energy = np.array(series.index)
        dos = np.array(series)
        if unit == "ev":
            energy, dos = dos_ry2ev(energy, series)
        return energy, dos

    def get_total_nos(self, spin=None, unit="ry"):
        unit = unit.lower()
        if unit not in UNITS:
            raise ValueError(f"Invalid unit: {unit}")

        df = self.get_tnos(spin, drop=True)
        series = df.groupby("E")["Total"].sum()
        energy = np.array(series.index)
        dos = np.array(series)
        if unit == "ev":
            energy, dos = dos_ry2ev(energy, dos)
        return energy, dos

    def get_partial_dos(
        self, sublatt=None, atom=None, spin=None, orbital="Total", unit="ry"
    ):
        unit = unit.lower()
        if unit not in UNITS:
            raise ValueError(f"Invalid unit: {unit}")

        df = self.get_pdos(sublatt, atom, spin, drop=True)
        series = df.groupby("E")[orbital].sum()
        energy = np.array(series.index)
        dos = np.array(series)
        if unit == "ev":
            energy, dos = dos_ry2ev(energy, dos)
        return energy, dos

    def get_total_slice(self, ita, it, s):
        dos = self.get_pdos()
        idx_ita = dos.index.unique("Sublatt")
        idx_it = dos.index.unique("Atom")
        idx_is = dos.index.unique("Spin")
        data = dos.loc[idx_ita[ita], idx_it[it], idx_is[s]]
        return data["Total"].array

    def dos_array(self, unit="ry"):
        """Return a 5D array of the total DOS with shape (ns, nit, nita, lmax, nzd)."""
        unit = unit.lower()
        if unit not in UNITS:
            raise ValueError(f"Invalid unit: {unit}")

        dos = self.get_pdos()
        idx_it = dos.index.unique("Sublatt")
        idx_ita = dos.index.unique("Atom")
        idx_ns = dos.index.unique("Spin")
        nita, nit, ns = len(idx_ita), len(idx_it), len(idx_ns)
        ncols = nita * nit * ns
        n = len(dos) // ncols
        columns = dos.columns[2:]  # Skip Energy and Total columns
        lmax = len(columns)
        shape = (ns, nit, nita, lmax, n)
        tdos = np.zeros(shape)
        energy = np.zeros(n)
        for i, ita in enumerate(idx_ita):
            for j, it in enumerate(idx_it):
                for k, s in enumerate(idx_ns):
                    data = dos.loc[it, ita, s]
                    energy[:] = data["E"].array
                    for c, col in enumerate(columns):
                        tdos[k, j, i, c] = data[col].array

        if unit == "ev":
            energy, dos = dos_ry2ev(energy, dos)
        return energy, tdos

    def total_dos_array(self, unit="ry"):
        unit = unit.lower()
        if unit not in UNITS:
            raise ValueError(f"Invalid unit: {unit}")

        dos = self.get_tdos()
        energy = np.array(dos["E"].array)
        dos = np.array(dos["Total"].array)
        if unit == "ev":
            energy, dos = dos_ry2ev(energy, dos)
        return energy, dos
