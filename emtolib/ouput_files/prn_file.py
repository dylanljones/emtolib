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

    def __init__(self, path):
        super().__init__(path)
        self.data = ""
        self.load()

    def loads(self, data: str) -> None:
        self.data = data

    def read_hopfield(self, unit="ev/aa^2"):
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
