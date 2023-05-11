# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

import re
import numpy as np
from datetime import datetime
from ..common import EmtoFile


RE_KEYVAL = re.compile(r"([a-zA-Z\(\)0-9]+).*?=.*?([a-zA-Z0-9_\-.]+)")

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

LAT2IBZ.update(dict((v, k) for k, v in LAT2IBZ.items()))

TEMPLATE = """\
BMDL      HP......={hp}                               {date}
JOBNAM...={jobname:10} MSGL.= {msgl:2d} MODE...={mode} NPRN.= {nprn:2d}
DIR001={dir001}
DIR006={dir006}
Madelung potential for {jobname}
NL.....= {nl:2d}
LAMDA....=    {lamda:6.2f} AMAX....=    {amax:6.2f} BMAX....=    {bmax:6.2f}
NQ....={nq:3d} LAT...={lat:2d} IPRIM.= {iprim} NQR2..= {nqr2}
A........= {a:9.7f} B.......= {b:9.7f} C.......= {c:9.7f}
"""


def lat_to_ibz(lat):
    """Returns the Bravais lattice IBZ code based on the input string code

    Parameters
    ----------
    lat : str
        The lattice name, for example 'bcc'.

    Returns
    -------
    ibz : int
        The Bravais lattice IBZ code, for example 3.
    """
    ltoi = {
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

    if lat not in ltoi.keys():
        raise ValueError(f"Unknown lattice {lat}!")

    return ltoi[lat]


class EmtoBmdlFile(EmtoFile):
    def __init__(self, path):
        now = datetime.now()
        super().__init__(path)
        self.jobname = None
        self.lat = None  # The lattice type (name)

        self.hp = "N"
        self.date = now.strftime("%d %b %y")
        self.msgl = 1
        self.mode = "B"  # B=Bulk, S=Surface
        self.nprn = 0
        self.nl = 7  # Number of orbitals in the Madelung matrix. (1-7)
        self.lamda = 2.5  # Constants in the Madelung matrix Ewald sum
        self.amax = 4.5
        self.bmax = 4.5
        self.nq = 1  # The number of atomic sites in the unit cell
        self.iprim = 1  # 1 if default primitive lattice vectors,  0 if user defined
        self.nqr2 = 0

        self.latconst = np.zeros(3, np.float64)
        self.latvectors = np.zeros(3, np.float64)
        self.latbasis = np.zeros(3, np.float64)

        self.dir001 = "bmdl/"
        self.dir006 = ""

    def loads(self, data: str):
        lines = data.splitlines(keepends=False)

        params = dict()
        # Line 1
        line = lines.pop(0)
        parts = line.split()
        params.update(dict(RE_KEYVAL.findall(parts[1])))
        params["date"] = " ".join(parts[2:])
        # Line 2
        params.update(dict(RE_KEYVAL.findall(lines.pop(0))))
        # Line 3+4
        self.dir001 = lines.pop(0).split("=")[1]
        self.dir006 = lines.pop(0).split("=")[1]
        lines.pop(0)  # Skip Line 5
        # Other lines
        params.update(dict(RE_KEYVAL.findall(" ".join(lines[:4]))))
        lines = lines[4:]
        self.hp = params["HP"]
        self.date = params["date"]
        self.jobname = params["JOBNAM"]
        self.msgl = int(params["MSGL"])
        self.nprn = int(params["NPRN"])
        self.nl = int(params["NL"])
        self.lamda = float(params["LAMDA"])
        self.amax = float(params["AMAX"])
        self.bmax = float(params["BMAX"])
        self.nq = int(params["NQ"])
        self.lat = LAT2IBZ[int(params["LAT"])]
        self.iprim = int(params["IPRIM"])
        self.nqr2 = int(params["NQR2"])

        self.latconst = np.array([float(params[k]) for k in ["A", "B", "C"]])
        if self.iprim:
            p = dict(RE_KEYVAL.findall(lines.pop(0)))
            self.latvectors = np.array([float(p[k]) for k in ["ALFA", "BETA", "GAMMA"]])
        else:
            self.latvectors = np.zeros((3, 3), np.float64)
            p = dict(RE_KEYVAL.findall(lines.pop(0)))
            for i in range(3):
                self.latvectors[i] = [float(p[k]) for k in ["BSX", "BSY", "BSZ"]]

        self.latbasis = np.zeros((self.nq, 3), np.float64)
        for i in range(self.nq):
            p = dict(RE_KEYVAL.findall(lines.pop(0)))
            self.latbasis[i, 0] = float(p[f"QX({i+1})"])
            self.latbasis[i, 1] = float(p[f"QY({i+1})"])
            self.latbasis[i, 2] = float(p[f"QZ({i+1})"])

    def dumps(self) -> str:
        params = {
            "hp": self.hp,
            "date": self.date,
            "jobname": self.jobname,
            "msgl": self.msgl,
            "mode": self.mode,
            "nprn": self.nprn,
            "dir001": self.dir001,
            "dir006": self.dir006,
            "nl": self.nl,
            "lamda": self.lamda,
            "amax": self.amax,
            "bmax": self.bmax,
            "nq": self.nq,
            "lat": LAT2IBZ[self.lat],
            "iprim": self.iprim,
            "nqr2": self.nqr2,
            "a": self.latconst[0],
            "b": self.latconst[1],
            "c": self.latconst[2],
        }
        data = TEMPLATE.format(**params)

        if self.iprim:
            line = "ALFA.....={0:10.6f} BETA....={1:10.6f} GAMMA...={2:10.6f}\n"
            data += line.format(*self.latvectors)
        else:
            line = "BSX......={0:10.7f} BSY......={1:10.7f} BSZ......={2:10.7f}\n"
            for i in range(3):
                data += line.format(*self.latvectors[i])

        line = "QX({0})....={1:10.7f} QY({0})...={2:10.7f} QZ({0})...={3:10.7f}\n"
        for i in range(self.nq):
            data += line.format(i + 1, *self.latbasis[i])
        return data
