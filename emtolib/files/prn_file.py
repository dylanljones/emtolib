# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-07

import re
import numpy as np
from ..common import EmtoFile, parse_params

RE_ATOM = re.compile("^Atom:(?P<atom>.*)")
RE_MAG = re.compile(r"^\s?Magn\. mom\. =\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
RE_MAG_ITER = re.compile(r"^\s?Magn\. mom\. =\s+(-?\d+\.\d+)$")
RE_DOS_EF = re.compile(r"DOS\(EF\)\s=\s+(?P<value>.*)")
RE_SECTION = re.compile(r"^ (?P<key>[A-Z]+):")
RE_ITER_LINE = re.compile(
    r"^\s*KGRN:\s+Iteration no\.(?P<iter>.*) Etot = (?P<etot>.*) erren = (?P<erren>.*)$"
)
RE_SWS_ALAT = re.compile(
    r"^\s*SWS =\s+(?P<sws>.*)\s+Alat =\s+(?P<alat_bohr>.*)\s+Bohr\s+(?P<alat>.*)\s+AA"
)
RE_ETOT = re.compile(r"^\s*Total energy\s*(?P<value>.*)$")
RE_ETOT_OKA = re.compile(r"^\s*Total energy\+OKA\s*(?P<value>.*)$")
RE_ETOT_EWALD = re.compile(r"^\s*Total\+Ewald\s*(?P<value>.*)$")


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

    text = prn.text
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


def parse_atom_panel(lines, istart):
    """Parse the atom panel from the PRN output file.

    Parameters
    ----------
    lines : list[str]
        The lines of the PRN file.
    istart : int
        The index of the first line of the atom panel.

    Returns
    -------
    data : dict
    """
    data = dict()

    i = istart
    # read header of panel:
    # Exch: ...
    while lines[i] != "":
        line = lines[i].strip()
        i += 1
        if line.startswith("Exch:"):
            exch = line.split()[0].split(":")[1]
            data["exch"] = exch
        elif line.startswith("Indices:"):
            indices = parse_params(line.split(" ", maxsplit=1)[1])
            data["it"] = int(indices["IT"])
            data["ita"] = int(indices["ITA"])

    i += 1
    # Read panel line
    #  Panel	    =  1   s -1/2 ...
    line = lines[i].strip()
    name, vals = line.split("=")
    assert name.strip() == "Panel"
    vals = [x.strip() for x in vals.split("  ") if x.strip()]
    data["panel"] = int(vals.pop(0))
    data["columns"] = vals

    i += 2
    # read first part of panel
    while lines[i] != "":
        line = lines[i]
        key, line = line.split("=")
        data[key.strip()] = [float(x) for x in line.split()]
        i += 1
    i += 1

    # Read second part of panel
    while lines[i] != "":
        line = lines[i]
        key, line = line.split("=")
        data[key.strip()] = [float(x) for x in line.split()]
        i += 1

    return data


class PrnFile(EmtoFile):

    extension = ".prn"
    autoload = True
    missing_ok = False

    def __init__(self, path):
        super().__init__(path)
        self.text = ""

    def loads(self, data: str) -> None:
        self.text = data

    @property
    def converged(self):
        return "Converged" in self.text

    def search_line(self, text, ignore_case=False):
        lines = list()
        for line in self.text.splitlines(keepends=False):
            if ignore_case:
                if text.lower() in line.lower():
                    lines.append(line)
            else:
                if text in line:
                    lines.append(line)
        return lines

    def grep(self, text, ignore_case=False):
        return "\n".join(self.search_line(text, ignore_case))

    def search(self, pattern, ignore_case=False):
        if ignore_case:
            return re.search(pattern, self.text, flags=re.IGNORECASE)
        else:
            return re.search(pattern, self.text)

    def findall(self, pattern, ignore_case=False):
        if ignore_case:
            return re.findall(pattern, self.text, flags=re.IGNORECASE)
        else:
            return re.findall(pattern, self.text)

    def iter_sections(self):
        """Iterate over the sections of the form ' KEY:' in the file.

        Yields
        ------
        key : str
            The key of the section (name of subroutine that wrote the section).
            The keys are *not* unique!
        lines : list of str
            The lines of the section.
        """
        lines = self.text.splitlines(keepends=False)
        i = 0
        start = -1
        last_key = ""
        while i < len(lines):
            line = lines[i]
            match = RE_SECTION.match(line)
            if match:
                key = match.group("key")
                if start >= 0:
                    end = i - 1
                    sec_lines = [line for line in lines[start:end] if line]
                    yield last_key, sec_lines
                last_key = key
                start = i
            i += 1
        sec_lines = [line for line in lines[start:] if line]
        yield last_key, sec_lines

    def get_sections(self, key):
        """Get all sections with the given key.

        Parameters
        ----------
        key : str
            The key of the sections (name of subroutine that wrote the section).

        Yields
        ------
        lines : list of str
            The lines of the sections with the given key.
        """
        for sec_key, sec_lines in self.iter_sections():
            if sec_key.lower() == key.lower():
                yield sec_lines

    def get_iterations(self):
        iterations = list()
        for line in self.search_line("Iteration"):
            match = RE_ITER_LINE.match(line)
            if match:
                it = {
                    "iter": int(match.group("iter")),
                    "etot": float(match.group("etot")),
                    "erren": float(match.group("erren")),
                }
                iterations.append(it)
        return iterations

    def get_lattice_constants(self):
        for line in self.search_line("SWS"):
            match = RE_SWS_ALAT.match(line)
            if match:
                sws = float(match.group("sws"))
                alat_bohr = float(match.group("alat_bohr"))
                alat = float(match.group("alat"))
                return sws, alat_bohr, alat
        return 0.0, 0.0, 0.0

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

    def get_magnetic_moment(self):
        pre = True
        current_atom = ""
        lines = self.text.splitlines(keepends=False)
        mag_pre, mag_post, mag_iter = list(), list(), list()
        for line in lines:
            line = line.strip()
            match = RE_ATOM.match(line)
            if match:
                current_atom = match["atom"]

            match = RE_MAG.match(line)
            if match:
                m1, m2 = float(match.group(1)), float(match.group(2))
                if pre:
                    mag_pre.append((current_atom, m1, m2))
                else:
                    mag_post.append((current_atom, m1, m2))
            match = RE_MAG_ITER.match(line)
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
        match = RE_DOS_EF.search(self.text)
        if match:
            return float(match.group("value"))
        return None

    def get_atom_panels(self):
        """Get the last atom output panels in the file."""
        panels = dict()
        regex = re.compile("^ Atom:(?P<atom>.*) S = (?P<s>.*) SWS = (?P<sws>.*)")
        lines = self.text.splitlines(keepends=False)
        for i, line in enumerate(lines):
            match = regex.match(line)
            if match:
                atom = match.group("atom").strip()
                s = float(match.group("s"))
                sws = float(match.group("sws"))
                data = parse_atom_panel(lines, i)
                data["sws"] = sws
                data["s"] = s
                panels[atom] = data
        return panels

    def get_etot(self):
        """Get the total energy from the last iteration."""
        iterations = self.get_iterations()
        return iterations[-1]["etot"]

    def get_total_energy(self):
        energies = dict()
        for line in self.search_line("Total"):
            match = RE_ETOT_EWALD.match(line)
            if match:
                try:
                    energies["etot+ewald"] = float(match.group("value"))
                except ValueError:
                    pass
                continue
            match = RE_ETOT_OKA.match(line)
            if match:
                try:
                    energies["etot+oka"] = float(match.group("value"))
                except ValueError:
                    pass
                continue
            match = RE_ETOT.match(line)
            if match:
                try:
                    energies["etot"] = float(match.group("value"))
                except ValueError:
                    pass
                continue
        if "etot" not in energies:
            last_iter = self.get_iterations()[-1]
            energies["etot"] = last_iter["etot"]
        return energies
