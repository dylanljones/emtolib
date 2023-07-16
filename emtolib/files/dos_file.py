# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import io
import re
import numpy as np
import pandas as pd
from typing import TextIO
from ..common import EmtoFile

RE_DOS_TOTAL = re.compile(r"Total DOS and NOS and partial \(IT\) (?P<spin>.*)$")
RE_DOS_SUBLAT = re.compile("Sublattice (?P<idx>.*) Atom (?P<atom>.*) spin (?P<spin>.*)")


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
            spin = match.group("spin").removeprefix("DOS")

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


class DosFile(EmtoFile):

    extensions = [".dos"]

    def __init__(self, *path):
        super().__init__(*path)
        self.tdos, self.pdos, self.tnos = load_dos(self.path)

    def get_tdos(self, spin=None, drop=True):
        return _index_dataframe(self.tdos, [spin], drop)

    def get_pdos(self, sublatt=None, atom=None, spin=None, drop=True):
        return _index_dataframe(self.pdos, [sublatt, atom, spin], drop)

    def get_tnos(self, sublatt=None, atom=None, spin=None, drop=True):
        return _index_dataframe(self.tnos, [sublatt, atom, spin], drop)

    def get_total_dos(self, spin=None):
        df = self.get_tdos(spin, drop=True)
        series = df.groupby("E")["Total"].sum()
        energy = np.array(series.index)
        dos = np.array(series)
        return energy, dos

    def get_total_nos(self, spin=None):
        df = self.get_tnos(spin, drop=True)
        series = df.groupby("E")["Total"].sum()
        energy = np.array(series.index)
        dos = np.array(series)
        return energy, dos

    def get_partial_dos(self, sublatt=None, atom=None, spin=None):
        df = self.get_pdos(sublatt, atom, spin, drop=True)
        series = df.groupby("E")["Total"].sum()
        energy = np.array(series.index)
        dos = np.array(series)
        return energy, dos
