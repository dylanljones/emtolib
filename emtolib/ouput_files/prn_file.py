# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import re
from ..common import EmtoFile


class EmtoPrnFile(EmtoFile):

    RE_ATOM = re.compile("^Atom:(?P<atom>.*)")

    RE_HOPFIELD = re.compile(r"^Hopfield\s+=\s+(?P<field1>.*),\s+(?P<field2>.*)")
    RE_FIELD = re.compile(r"(?P<value>.*)\((?P<unit>.*)\)")

    RE_MAG = re.compile(r"^\s?Magn\. mom\. =\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
    RE_MAG_ITER = re.compile(r"^\s?Magn\. mom\. =\s+(-?\d+\.\d+)$")

    def __init__(self, path):
        super().__init__(path)
        self.data = ""
        self.load()

    def loads(self, data: str) -> None:
        self.data = data

    @property
    def converged(self):
        return "Converged" in self.data

    def get_hopfield(self, unit="ev/aa^2"):
        # Prepare and check unit
        unit = unit.lower()
        if unit.startswith("ev"):
            unit = "ev/aa^2"
        elif unit.startswith("ry"):
            unit = "ry/bohr^2"
        if unit not in ("ry/bohr^2", "ev/aa^2"):
            raise ValueError(f"Unit must be 'Ry/Bohr^2' or 'eV/AA^2', not '{unit}'")

        fields = dict()

        current_atom = ""
        lines = self.data.splitlines(keepends=False)
        for line in lines:
            line = line.strip()
            match = self.RE_ATOM.match(line)
            if match:
                current_atom = match["atom"]

            match = self.RE_HOPFIELD.match(line)
            if match:
                s1 = match["field1"]
                s2 = match["field2"]
                match = self.RE_FIELD.match(s1)
                field1 = float(match["value"])
                unit1 = match["unit"].lower()
                if unit == unit1:
                    fields[current_atom] = field1
                    continue
                match = self.RE_FIELD.match(s2)
                field2 = float(match["value"])
                unit2 = match["unit"].lower()
                assert unit == unit2
                fields[current_atom] = field2
        return fields

    def read_hopfield(self, unit="ev/aa^2"):
        """Backwards compatibility"""
        return self.get_hopfield(unit)

    def read_magnetic_moment(self):
        pre = True
        current_atom = ""
        lines = self.data.splitlines(keepends=False)
        mag_pre, mag_post, mag_iter = list(), list(), list()
        for line in lines:
            line = line.strip()
            match = self.RE_ATOM.match(line)
            if match:
                current_atom = match["atom"]

            match = self.RE_MAG.match(line)
            if match:
                m1, m2 = float(match.group(1)), float(match.group(2))
                if pre:
                    mag_pre.append((current_atom, m1, m2))
                else:
                    mag_post.append((current_atom, m1, m2))
            match = self.RE_MAG_ITER.match(line)
            if match:
                pre = False
                mag_iter.append((current_atom, float(match.group(1))))

        return mag_pre, mag_post, mag_iter

    def get_total_magnetic_moment(self):
        mag_pre, mag_post, mag_iter = self.read_magnetic_moment()
        mtot = mag_post[0][-1]
        for x in mag_post[1:]:
            assert x[-1] == mtot
        return mtot

    def get_atomic_magnetic_moment(self, pre=False):
        mag_pre, mag_post, mag_iter = self.read_magnetic_moment()
        moments = mag_pre if pre else mag_post
        return [(at, m) for at, m, _ in moments]
