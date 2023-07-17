# coding: utf-8
#
# This code is part of NbTa_Superconductor.
#
# Copyright (c) 2022, Dylan Jones

import re
import numpy as np
from ..common import EmtoFile


def extract_hopfields(prn, unit="ev/aa^2"):
    """Extracts the Hopfield parameters in the given unit from the given PRN file.

    Parameters
    ----------
    prn : PrnFile
        The PRN file to extract the Hopfield parameters from.
    unit : str, optional
        The unit of the Hopfield parameters. Can be 'Ry/Bohr^2' or 'eV/AA^2'.

    Returns
    -------
    hopfields : list[tuple[str, list[float]]]
        A list of tuples. Each tuple contains the atom symbol and a list of the
        Hopfield parameters for that atom and for each spin
    """
    re_atom = re.compile("^Atom:(?P<atom>.*)")
    re_hopfield = re.compile(r"^Hopfield\s+=\s+(?P<field1>.*),\s+(?P<field2>.*)")

    # Prepare and check unit
    unit = unit.lower()
    if unit.startswith("ev"):
        unit = "ev/aa^2"
    elif unit.startswith("ry"):
        unit = "ry/bohr^2"
    if unit not in ("ry/bohr^2", "ev/aa^2"):
        raise ValueError(f"Unit must be 'Ry/Bohr^2' or 'eV/AA^2', not '{unit}'")

    text = prn.data
    hopfields = list()
    current_atom = ""
    ihf = 0 if unit == "ry/bohr^2" else 1

    hopfield = list()

    # Hopfields are stored in blocks:
    # Atom:<Symbol>
    # ...
    # Hopfield = 0.013988 (Ry / Bohr ^ 2), 0.679618 (eV / AA ^ 2)
    # ...
    for line in text.splitlines(keepends=False):
        line = line.strip()

        # Check for new Atom line
        match = re_atom.match(line)
        if match:
            if hopfield:
                # New block: We already have seen hopfields, this means a new block
                # begins here. Store the hopfields of the previous block.
                hopfields.append((current_atom, hopfield))
                hopfield = list()
            # Update current atom
            current_atom = match["atom"]

        # Check line for Hopfield
        match = re_hopfield.match(line)
        if match:
            # Extract the hopfield values
            etas = [match["field1"], match["field2"]]
            # Choose the correct unit and convert to float
            eta = float(etas[ihf].split("(")[0])
            # Add to list of hopfield values of *current* block
            hopfield.append(eta)

    if hopfield:
        # Store the last block
        hopfields.append((current_atom, hopfield))

    return hopfields


class PrnFile(EmtoFile):

    extension = ".prn"

    RE_ATOM = re.compile("^Atom:(?P<atom>.*)")

    RE_HOPFIELD = re.compile(r"^Hopfield\s+=\s+(?P<field1>.*),\s+(?P<field2>.*)")
    RE_FIELD = re.compile(r"(?P<value>.*)\((?P<unit>.*)\)")

    RE_MAG = re.compile(r"^\s?Magn\. mom\. =\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
    RE_MAG_ITER = re.compile(r"^\s?Magn\. mom\. =\s+(-?\d+\.\d+)$")

    RE_DOS_EF = re.compile(r"DOS\(EF\)\s=\s+(?P<value>.*)")

    def __init__(self, path):
        super().__init__(path)
        self.data = ""
        self.load(missing_ok=True)

    def loads(self, data: str) -> None:
        self.data = data

    @property
    def converged(self):
        return "Converged" in self.data

    def search_line(self, text, ignore_case=False):
        lines = list()
        for line in self.data.splitlines(keepends=False):
            if ignore_case:
                if text.lower() in line.lower():
                    lines.append(line)
            else:
                if text in line:
                    lines.append(line)
        return lines

    def grep(self, text, ignore_case=False):
        return "\n".join(self.search_line(text, ignore_case))

    def search(self, pattern):
        return re.search(pattern, self.data)

    def findall(self, pattern):
        return re.findall(pattern, self.data)

    def extract_hopfields(self, unit="ev/aa^2"):
        return extract_hopfields(self, unit)

    def get_sublat_hopfields(self, dat, unit="ev/aa^2"):
        # Extract raw Hopfield parameters
        raw = self.extract_hopfields(unit)
        assert len(raw) == len(dat.atoms)
        assert tuple(x[0] for x in raw) == tuple(x.symbol for x in dat.atoms)

        # Construct sublattice indices
        sublattices = {i: list() for i in range(dat.nt)}
        for i, atom in enumerate(dat.atoms):
            sublattices[atom.it - 1].append(i)

        # Average over sublattices
        hopfields = list()
        for it in range(dat.nt):
            # Get the Hopfield parameters for each atom in the sublattice
            etas = np.array([raw[i][1] for i in sublattices[it]])
            # Get the concentrations of each atom in the sublattice
            concs = [dat.atoms[i].conc for i in sublattices[it]]
            # Average over the sublattice
            eta = np.average(etas, weights=concs, axis=0)
            hopfields.append(eta)

        return np.array(hopfields)

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

    def get_magnetic_moment(self):
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
        mag_pre, mag_post, mag_iter = self.get_magnetic_moment()
        mtot = mag_post[0][-1]
        for x in mag_post[1:]:
            assert x[-1] == mtot
        return mtot

    def get_atomic_magnetic_moment(self, pre=False):
        mag_pre, mag_post, mag_iter = self.get_magnetic_moment()
        moments = mag_pre if pre else mag_post
        return [(at, m) for at, m, _ in moments]

    def get_dos_ef(self):
        match = self.RE_DOS_EF.search(self.data)
        if match:
            return float(match.group("value"))
        return None
