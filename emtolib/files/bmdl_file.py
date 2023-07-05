# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

from datetime import datetime
import numpy as np
from ..common import EmtoFile, parse_params
from ..ftmplt import Template

TEMPLATE = """\
BMDL      HP......={hp} {header:>39}
JOBNAM...={jobnam:10} MSGL.= {msgl:2d} NPRN.= {nprn:2d}
DIR001={dir001}
DIR006={dir006}
{comment}
NL.....= {nl:d}
LAMDA....=    {lamda:6.2f} AMAX....=    {amax:6.2f} BMAX....=    {bmax:6.2f}
{cell}
"""

# Lattice name to IBZ code
LAT2IBZ = {
    "sc": 1,
    "fcc": 2,
    "bcc": 3,
    "hcp": 4,
    "st": 5,
    "bct": 6,
    "trig": 7,
    "so": 8,
    "baco": 9,
    "bco": 10,
    "fco": 11,
    "sm": 12,
    "bacm": 13,
    "bcm": 13,
    "stric": 14,
    "B2": 1,
    "L12": 1,
}
# IBZ code to lattice name
IBZ2LAT = dict((v, k) for k, v in LAT2IBZ.items())


def lat_to_ibz(lat):
    """Returns the Bravais lattice IBZ code based on the input string code"""
    if lat not in LAT2IBZ.keys():
        valid = list(LAT2IBZ.keys())
        raise ValueError(f"Unknown lattice {lat}! Valid lattices: {valid}")
    return LAT2IBZ[lat]


def ibz_to_lat(ibz):
    """Returns the Bravais lattice name based on the input IBZ code"""
    if ibz not in IBZ2LAT.keys():
        valid = list(IBZ2LAT.keys())
        raise ValueError(f"Unknown IBZ code {ibz}! IBZ codes: {valid}")
    return IBZ2LAT[ibz]


def _parse_line(line, dtype=float):
    values, info = list(), ""
    parts = line.split("=")

    for part in parts[1:-1]:
        val = part.split()[0]
        values.append(dtype(val))
    items = parts[-1].split()
    values.append(dtype(items[0]))
    if len(items) > 1:
        info = items[-1]
    return values, info


def parse_cell(cell, raise_errors=True):
    lines = cell.splitlines()
    params = parse_params(lines.pop(0))
    constants, _ = _parse_line(lines.pop(0))
    basis_vecs, basis_vecs_info = list(), list()
    latt_vecs, latt_vecs_info = list(), list()
    angles = list()
    for line in lines:
        if line.upper().startswith("Q"):
            # Positions in the conventional unit cell
            vals, info = _parse_line(line)
            basis_vecs.append(vals)
            basis_vecs_info.append(info)
        elif line.upper().startswith("B"):
            # Lattice vectors in units of A
            vals, info = _parse_line(line)
            basis_vecs.append(vals)
            basis_vecs_info.append(info)
        elif line.upper().startswith("AL"):
            angles, _ = _parse_line(line)
        else:
            if raise_errors:
                raise ValueError(f"Unknown line: {line}")

    params["constants"] = constants
    params["angles"] = angles
    params["latt_vecs"] = latt_vecs
    params["latt_vecs_info"] = latt_vecs_info
    params["basis_vecs"] = basis_vecs
    params["basis_vecs_info"] = basis_vecs_info

    return params


class EmtoBmdlFile(EmtoFile):

    extension = ".dat"
    template = Template(TEMPLATE, ignore_case=True)

    def __init__(self, path, **kwargs):
        super().__init__(path)

        self.jobnam = "kgrn"
        self.hp = "N"
        self.header = ""  # first line after KGRN (usually date in the format %d %b %y)
        self.msgl = 1  # Determines what/how much will be printed on screen
        self.nprn = 0  # Determines what/how much will be printed in the output file
        self.dir001 = "mdl/"  # Directory, where the slope matrix *.tfh will be stored
        self.dir006 = ""  # Directory, where the output file *.prn will be stored.

        # Line after path variables, eg 'Self-consistent KKR calculation for {jobnam}'
        self.comment = "Madelung potential for {jobnam}"

        self.nl = 0  # Number of orbitals (possible choices: 1, 2, 3, 4 or 5).
        self.lamda = 2.5  # Constants in the Madelung matrix Ewald sum
        self.amax = 4.5  # The maximum radius of the atomic spheres
        self.bmax = 4.5  # The maximum radius of the muffin-tin spheres

        self.nq = 0
        self.nq3 = None
        self.nghbp = None
        self.lat = "sc"
        self.iprim = 0
        self.nqr2 = 0

        self.constants = list()
        self.angles = None
        self.lat_vecs = None
        self.lat_vecs_info = None
        self.basis_vecs = list()
        self.basis_vecs_info = list()

        self.load(missing_ok=True)
        if kwargs:
            self.update(kwargs)

    def get_lattice_vectors(self, scaled=True):
        """Returns the lattice vectors in units of A

        If the BMDL file contains lattice vectors, those are returned.
        Otherwise, the lattice vectors are computed from the angles:
        The first vector is along the x-axis
        and the second vector is in the xy-plane. The third vector is computed
        using all three angles. The angles are the opposite angles of the vectors,
        so alpha is the angles between a_2 and a_3, beta is the angle between
        a_1 and a_3 and gamma is the angles beteen a_1 and a_2.
              a_3
               |
               | α
             β +------>a_2
              ╱ γ
             ╱
           a_1

        Parameters
        ----------
        scaled : bool
            If True, the lattice vectors are scaled by the lattice constants.

        Returns
        -------
        (N, N) np.ndarray
            Lattice vectors in units of A. The vectors are the rows of the matrix.
        """
        if self.lat_vecs is not None:
            vectors = self.lat_vecs
        else:
            alpha, beta, gamma = np.deg2rad(self.angles)
            v1 = np.array([1, 0, 0])
            v2 = np.array([np.cos(gamma), np.sin(gamma), 0])
            x = np.cos(beta)
            y = np.cos(alpha)
            z = np.sin(beta) * np.sin(alpha)
            v3 = np.array([x, y, z])
            vectors = np.array([v1, v2, v3])
        if scaled:
            vectors = (vectors.T * self.constants).T
        return vectors

    def transform(self, positions):
        r"""Transforms positions from fractional to cartesian coordinates

        The given positions .math:`x_i` are transformed to cartesian coordinates
        using the lattice vectors .math:`A = \sum_{ij} a_{ij}` as
        .. math::
            x'_i = \sum_j x_j a_{ij}

        Parameters
        ----------
        positions : (..., N) array_like
            Positions in fractional coordinates

        Returns
        -------
        (..., N) ndarray
            Positions in cartesian coordinates
        """
        positions = np.asarray(positions)
        vectors = self.get_lattice_vectors(scaled=True)
        trans = np.dot(positions, vectors.T)
        assert trans.shape == positions.shape
        return trans

    def itransform(self, positions):
        """Transforms positions from cartesian to fractional coordinates

        Parameters
        ----------
        positions : (..., N) array_like
            Positions in fractional coordinates

        Returns
        -------
        (..., N) ndarray
            Positions in cartesian coordinates
        """
        positions = np.asarray(positions)
        vectors = self.get_lattice_vectors(scaled=True)
        inv_vectors = np.linalg.inv(vectors)
        trans = np.dot(positions, inv_vectors.T)
        assert trans.shape == positions.shape
        return trans

    def get_site_ratio(self, tol=1e-7):
        """Computes the ration of a site in the unit cell."""
        positions = np.asarray(self.basis_vecs)
        # get fractional coordinates (inversed basis transformation)
        coord = self.itransform(positions)
        coord = coord % 1  # ensure the coordinates are between 0 and 1
        # fractional coordinates on the boundary of the cell, i.e. 0 or 1
        bound = np.isclose(coord, 0, atol=tol) | np.isclose(coord, 1, atol=tol)
        # number of coordinates on edges:
        # 0 (interior): site in 1 cell
        # 1 (face):     site in 2 cells
        # 2 (edge):     site in 4 cells
        # 3 (corner):   site in 8 cells
        ncell = 2 ** np.count_nonzero(bound, axis=-1)
        # Check output shape
        assert ncell.shape == positions.shape[:-1]
        return 1 / ncell

    def set_header(self, header="", date_frmt="%d %b %y"):
        self.header = (header + " " + datetime.now().strftime(date_frmt)).strip()

    def update(self, *args, **kwargs):
        data = dict(*args, **kwargs)
        for k, v in data.items():
            if not hasattr(self, k):
                raise KeyError(f"{k} is not a valid field of {self.__class__.__name__}")
            if isinstance(v, str):
                v = v.strip()
            self.__setattr__(k, v)

    def __getitem__(self, key):
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        return self.__getattribute__(key)

    def __setitem__(self, key, value):
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        self.__setattr__(key, value)

    def to_dict(self):
        data = {k: v for k, v in self.__dict__.items() if not k.startswith("_")}
        data.pop("cell", None)
        return data

    def check(self):
        if self.nq is not None:
            if len(self.basis_vecs) != self.nq:
                raise ValueError(
                    f"Number of basis vectors ({len(self.basis_vecs)}) "
                    f"does not match NQ ({self.nq})"
                )
        elif self.nq3 is not None:
            if len(self.basis_vecs) != self.nq3:
                raise ValueError(
                    f"Number of basis vectors ({len(self.basis_vecs)}) "
                    f"does not match NQ3 ({self.nq3})"
                )
        else:
            raise ValueError("Must specify either NQ or NQ3")

        if self.angles is not None and self.lat_vecs is not None:
            raise ValueError("Cannot specify both angles and lattice vectors")

    def loads(self, text: str) -> None:
        data = self.template.parse(text.replace(".d", ".e"))
        params = dict(data.copy())
        cell_params = parse_cell(params.pop("cell"), raise_errors=True)
        self.update(params)

        try:
            self.nq = int(cell_params.pop("NQ"))
            self.nq3 = None
        except KeyError:
            self.nq = None
            self.nq3 = int(cell_params.pop("NQ3"))
        self.lat = ibz_to_lat(int(cell_params.pop("LAT")))
        self.iprim = int(cell_params.pop("IPRIM"))
        try:
            self.nghbp = int(cell_params.pop("NHGBP"))
        except KeyError:
            self.nghbp = None
        self.nqr2 = int(cell_params.pop("NQR2"))

        self.constants = np.array(cell_params.pop("constants"))
        angles = cell_params.pop("angles")
        if angles:
            self.angles = np.array(angles)
        else:
            self.angles = None
        lat_vecs = cell_params.pop("latt_vecs")
        lat_vecs_info = cell_params.pop("latt_vecs_info")
        if lat_vecs:
            self.lat_vecs = np.array(lat_vecs)
            self.lat_vecs_info = lat_vecs_info
        else:
            self.lat_vecs = None
            self.lat_vecs_info = None
        self.basis_vecs = np.array(cell_params.pop("basis_vecs"))
        self.basis_vecs_info = cell_params.pop("basis_vecs_info")

        self.check()

    def dumps(self) -> str:
        raise NotImplementedError
