# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-21

import re
import logging
from datetime import datetime
from typing import Union, List, Set, Dict
from pathlib import Path
import numpy as np
from ..common import EmtoFile, parse_params, elements, dict_diff
from ..errors import KGRNError, KGRNReadError, KGRNWriteError
from ..ftmplt import Template
from ..config import update_emto_paths

logger = logging.getLogger(__name__)

RE_BAND_SECTION = re.compile(r"Band: (.*?) lines")

TEMPLATE_KGRN = """\
KGRN {header:>55}
JOBNAM={jobnam}
STRT..=  {strt} MSGL.={msgl:3d} EXPAN.= {expan} FCD..=  {fcd} FUNC..= {func}
FOR001={for001}
FOR001={for001_2}
DIR002={dir002}
DIR003={dir003}
FOR004={for004}
DIR006={dir006}
DIR009={dir009}
DIR010={dir010}
DIR011={dir011}
{comment}
Band: 10 lines
NITER.={niter:3d} NLIN.={nlin:3d} NPRN.={nprn:3d} NCPA.={ncpa:3d} NT...={nt:3d} MNTA.={mnta:3d}
MODE..= {mode:2} FRC..=  {frc} DOS..=  {dos} OPS..=  {ops} AFM..=  {afm} CRT..=  {crt}
Lmaxh.={lmaxh:3d} Lmaxt={lmaxt:3d} NFI..={nfi:3d} FIXG.={fixg:3d} SHF..={shf:3d} SOFC.=  {sofc}
KMSH...= {kmsh} IBZ..={ibz:3d} NKX..={nkx:3d} NKY..={nky:3d} NKZ..={nkz:3d} FBZ..=  {fbz}
KMSH2..= {kmsh2} IBZ2.={ibz2:3d} NKX2.={nkx2:3d} NKY2.={nky2:3d} NKZ2.={nkz2:3d}
ZMSH...= {zmsh} NZ1..={nz1:3d} NZ2..={nz2:3d} NZ3..={nz3:3d} NRES.={nres:3d} NZD.={nzd:4d}
DEPTH..={depth:7.3f} IMAGZ.={imagz:7.5f} EPS...={eps:7.3f} ELIM..={elim:7.3f}
AMIX...={amix:7.5f} EFMIX.={efmix:7.5f} VMTZ..={vmtz:7.3f} MMOM..={mmom:7.3f}
TOLE...= {tole:7.1e} TOLEF.= {tolef:7.1e} TOLCPA= {tolcpa:7.1e} TFERMI={tfermi:7.3f} (K)
SWS......={sws:10.6f} NSWS.={nsws:3d} DSWS..={dsws:7.2f} ALPCPA={alpcpa:7.4f}
Setup: 2 + NQ*NS lines
EFGS...={efgs:7.3f} HX....={hx:7.3f} NX...={nx:3d} NZ0..={nz0:3d} STMP..= {stmp}
{atoms}
Atom:  4 lines + NT*NTA*6 lines
IEX...={iex:3d} NP..= {np:3d} NES..={nes:3d} NITER={dirac_niter:3d} IWAT.={iwat:3d} NPRNA={nprna:3d}
VMIX.....={vmix:10.6f} RWAT....={rwat:10.6f} RMAX....={rmax:10.6f}
DX.......={dx:10.6f} DR1.....={dr1:10.6f} TEST....=  {test:8.2E}
TESTE....=  {teste:8.2E} TESTY...=  {testy:8.2E} TESTV...=  {testv:8.2E}
{atomconf}
"""  # noqa

TEMPLATE_DMFT = """\
DMFT {header:>55}
JOBNAM={jobnam}
STRT..=  {strt} MSGL.={msgl:3d} EXPAN.= {expan} FCD..=  {fcd} FUNC..= {func}
FOR001={for001}
FOR001={for001_2}
DIR002={dir002}
DIR003={dir003}
FOR004={for004}
FOR007={for007}
DIR006={dir006}
DIR009={dir009}
DIR010={dir010}
DIR011={dir011}
{comment}
Band: 10 lines
NITER.={niter:3d} NLIN.={nlin:3d} NPRN.={nprn:3d} NCPA.={ncpa:3d} NT...={nt:3d} MNTA.={mnta:3d}
MODE..= {mode:2} FRC..=  {frc} DOS..=  {dos} OPS..=  {ops} AFM..=  {afm} CRT..=  {crt}
Lmaxh.={lmaxh:3d} Lmaxt={lmaxt:3d} NFI..={nfi:3d} FIXG.={fixg:3d} SHF..={shf:3d} SOFC.=  {sofc}
KMSH...= {kmsh} IBZ..={ibz:3d} NKX..={nkx:3d} NKY..={nky:3d} NKZ..={nkz:3d} FBZ..=  {fbz}
KMSH2..= {kmsh2} IBZ2.={ibz2:3d} NKX2.={nkx2:3d} NKY2.={nky2:3d} NKZ2.={nkz2:3d}
ZMSH...= {zmsh} NZ1..={nz1:3d} NZ2..={nz2:3d} NZ3..={nz3:3d} NRES.={nres:3d} NZD.={nzd:4d}
DEPTH..={depth:7.3f} IMAGZ.={imagz:7.5f} EPS...={eps:7.3f} ELIM..={elim:7.3f}
AMIX...={amix:7.5f} EFMIX.={efmix:7.5f} VMTZ..={vmtz:7.3f} MMOM..={mmom:7.3f}
TOLE...= {tole:7.1e} TOLEF.= {tolef:7.1e} TOLCPA= {tolcpa:7.1e} TFERMI={tfermi:7.3f} (K)
SWS......={sws:10.6f} NSWS.={nsws:3d} DSWS..={dsws:7.2f} ALPCPA={alpcpa:7.4f}
NOM...={nom:4d} NOMI.={nomi:4d} DC.={dc:3d} TTT...={ttt:7.3f} SMIX..={smix:7.3f}
SOLVER={solver}
Setup: 2 + NQ*NS lines
EFGS...={efgs:7.3f} HX....={hx:7.3f} NX...={nx:3d} NZ0..={nz0:3d} STMP..= {stmp}
{atoms}
Atom:  4 lines + NT*NTA*6 lines
IEX...={iex:3d} NP..= {np:3d} NES..={nes:3d} NITER={dirac_niter:3d} IWAT.={iwat:3d} NPRNA={nprna:3d}
VMIX.....={vmix:10.6f} RWAT....={rwat:10.6f} RMAX....={rmax:10.6f}
DX.......={dx:10.6f} DR1.....={dr1:10.6f} TEST....=  {test:8.2E}
TESTE....=  {teste:8.2E} TESTY...=  {testy:8.2E} TESTV...=  {testv:8.2E}
{atomconf}
"""  # noqa

TEMPLATE_DMFT_v2 = """\
DMFT {header:>55}
JOBNAM={jobnam}
STRT..=  {strt} MSGL.={msgl:3d} EXPAN.= {expan} FCD..=  {fcd} FUNC..= {func}
FOR001={for001}
FOR001={for001_2}
DIR002={dir002}
DIR003={dir003}
FOR004={for004}
FOR007={for007}
DIR006={dir006}
DIR009={dir009}
DIR010={dir010}
DIR011={dir011}
{comment}
Band: 10 lines
NITER.={niter:3d} NLIN.={nlin:3d} NPRN.={nprn:3d} NCPA.={ncpa:3d} NT...={nt:3d} MNTA.={mnta:3d}
MODE..= {mode:2} FRC..=  {frc} DOS..=  {dos} OPS..=  {ops} AFM..=  {afm} CRT..=  {crt}
Lmaxh.={lmaxh:3d} Lmaxt={lmaxt:3d} NFI..={nfi:3d} FIXG.={fixg:3d} SHF..={shf:3d} SOFC.=  {sofc}
KMSH...= {kmsh} IBZ..={ibz:3d} NKX..={nkx:3d} NKY..={nky:3d} NKZ..={nkz:3d} FBZ..=  {fbz}
KMSH2..= {kmsh2} IBZ2.={ibz2:3d} NKX2.={nkx2:3d} NKY2.={nky2:3d} NKZ2.={nkz2:3d} AVG..=  {avg}
ZMSH...= {zmsh} NZ1..={nz1:3d} NZ2..={nz2:3d} NZ3..={nz3:3d} NRES.={nres:3d} NZD.={nzd:4d}
DEPTH..={depth:7.3f} IMAGZ.={imagz:7.5f} EPS...={eps:7.3f} ELIM..={elim:7.3f}
AMIX...={amix:7.5f} EFMIX.={efmix:7.5f} VMTZ..={vmtz:7.3f} MMOM..={mmom:7.3f}
TOLE...= {tole:7.1e} TOLEF.= {tolef:7.1e} TOLCPA= {tolcpa:7.1e} TFERMI={tfermi:7.3f} (K)
SWS......={sws:10.6f} NSWS.={nsws:3d} DSWS..={dsws:7.2f} ALPCPA={alpcpa:7.4f}
NOM...={nom:4d} NOMI.={nomi:4d} DC.={dc:3d} TTT...={ttt:7.3f} SMIX..={smix:7.3f}
SOLVER={solver}
Setup: 2 + NQ*NS lines
EFGS...={efgs:7.3f} HX....={hx:7.3f} NX...={nx:3d} NZ0..={nz0:3d} STMP..= {stmp}
{atoms}
Atom:  4 lines + NT*NTA*6 lines
IEX...={iex:3d} NP..= {np:3d} NES..={nes:3d} NITER={dirac_niter:3d} IWAT.={iwat:3d} NPRNA={nprna:3d}
VMIX.....={vmix:10.6f} RWAT....={rwat:10.6f} RMAX....={rmax:10.6f}
DX.......={dx:10.6f} DR1.....={dr1:10.6f} TEST....=  {test:8.2E}
TESTE....=  {teste:8.2E} TESTY...=  {testy:8.2E} TESTV...=  {testv:8.2E}
{atomconf}
"""  # noqa


ATOM_COLUMNS = tuple("iq it ita nz conc sms sws wswst qtr splt fix".split())
ATCONF_TEMPLATE = "Iz= {iz:3d} Norb={norb:3d} Ion=  {ion:d} Config= {config}"
ATLINE_TEMPLATE = (
    "{symbol:2}    {iq:3d} {it:2d} {ita:2d}  {nz:2d}  {conc:5.3f}  "
    "{sms:5.3f}  {sws:5.3f}  {wswst:5.3f} {qtr:4.1f} {splt:4.1f}  {fix}"
)
ORBITALS = ["s", "p", "d", "f", "g", "h"]

DEFAULT_PARAMS = {
    "jobnam": "nb",
    "strt": "A",
    "msgl": 1,
    "expan": "S",
    "fcd": "Y",
    "func": "SCA",
    "for001": "~/EMTO/EMTO5.8/kstr/smx/bcc.tfh",
    "for001_2": "",
    "for004": "~/EMTO/EMTO5.8/bmdl/mdl/bcc.mdl",
    "dir002": "pot/",
    "dir003": "pot/",
    "dir006": "",
    "dir009": "pot/",
    "dir010": "chd/",
    "dir011": "/tmp/",
    "niter": 500,
    "nlin": 31,
    "nprn": 0,
    "ncpa": 30,
    "nt": 1,
    # "mnta": 2,
    "mode": "3D",
    "frc": "N",
    "dos": "D",
    "ops": "N",
    "afm": "P",
    "crt": "M",
    "lmaxh": 8,
    "lmaxt": 4,
    "nfi": 31,
    "fixg": 2,
    "shf": 0,
    "sofc": "N",
    "kmsh": "G",
    "ibz": 3,
    "nkx": 0,
    "nky": 75,
    "nkz": 0,
    "fbz": "N",
    "kmsh2": "G",
    "ibz2": 1,
    "nkx2": 4,
    "nky2": 0,
    "nkz2": 51,
    "zmsh": "C",
    "nz1": 32,
    "nz2": 8,
    "nz3": 8,
    "nres": 4,
    "nzd": 5000,
    "depth": 1.0,
    "imagz": 0.01,
    "eps": 0.2,
    "elim": -1.0,
    "amix": 0.1,
    "efmix": 1.0,
    "vmtz": 0.0,
    "mmom": 0.0,
    "tole": 1e-7,
    "tolef": 1e-7,
    "tolcpa": 1e-6,
    "tfermi": 500.0,
    "sws": 2.8388,
    "nsws": 1,
    "dsws": 0.05,
    "alpcpa": 0.602,
    # ----------------
    "efgs": 0.0,
    "hx": 0.3,
    "nx": 9,
    "nz0": 16,
    "stmp": "N",
    # ----------------
    "iex": 4,
    "np": 251,
    "nes": 15,
    "dirac_niter": 100,
    "iwat": 0,
    "nprna": 0,
    "vmix": 0.3,
    "rwat": 3.5,
    "rmax": 20.0,
    "dx": 0.03,
    "dr1": 0.002,
    "test": 1e-12,
    "teste": 1e-12,
    "testy": 1e-12,
    "testv": 1e-12,
}


def format_atom_line(params):
    line = ATLINE_TEMPLATE.format(**params)
    uj = params["u"] + params["j"]
    if uj:
        line += " " + " ".join(f"{x:4.2f}" for x in uj)
    return line.strip()


def format_atom_block(params):
    lines = list()
    lines.append(params["symbol"])
    lines.append(ATCONF_TEMPLATE.format(**params))
    lines.append("n     " + " ".join(f"{x:2d}" for x in params["n"]))
    lines.append("Kappa " + " ".join(f"{x:2d}" for x in params["kappa"]))
    lines.append("Occup " + " ".join(f"{x:2d}" for x in params["occup"]))
    lines.append("Valen " + " ".join(f"{x:2d}" for x in params["valen"]))
    return "\n".join(lines)


def parse_atom_line(line, columns=ATOM_COLUMNS):
    ncols = len(columns)
    parts = line.split()
    at = parts.pop(0)
    values = parts[:ncols]
    params = dict(zip(columns, values))
    params["iq"] = int(params["iq"])
    params["it"] = int(params["it"])
    params["ita"] = int(params["ita"])
    params["nz"] = int(params["nz"])
    params["conc"] = float(params["conc"])
    params["sms"] = float(params["sms"])
    params["sws"] = float(params["sws"])
    params["wswst"] = float(params["wswst"])
    params["qtr"] = float(params["qtr"])
    params["splt"] = float(params["splt"])
    if "fix" not in params:
        params["fix"] = ""

    uj = [float(x) for x in parts[ncols:]]
    nuj = len(uj)
    params["u"] = uj[: nuj // 2]
    params["j"] = uj[nuj // 2 :]
    return at, params


def parse_atom_block(lines):
    assert len(lines) == 6
    at = lines.pop(0).strip()
    at_params = {k.lower(): v for k, v in parse_params(lines.pop(0).strip()).items()}
    at_params["iz"] = int(at_params["iz"])
    at_params["norb"] = int(at_params["norb"])
    at_params["ion"] = int(at_params["ion"])
    for line in lines:
        parts = line.strip().split()
        key = parts.pop(0)
        values = [int(x) for x in parts]
        at_params[key.lower()] = values
    return at, at_params


def parse_atoms(atomstr, atomconfstr):
    atom_lines = atomstr.splitlines()
    header = atom_lines.pop(0).lower()
    columns = header.replace("(", "").replace(")", "").split()[1:]
    if tuple(columns) != ATOM_COLUMNS[: len(columns)]:
        raise KGRNReadError(f"Invalid atom header: {header}")
    atom_params = [parse_atom_line(line, columns) for line in atom_lines]

    lines = atomconfstr.splitlines()
    n_blocks = len(lines) // 6
    atom_blocks = [
        parse_atom_block(lines[i * 6 : (i + 1) * 6]) for i in range(n_blocks)
    ]

    counts = dict()
    atoms = list()
    for at, params in atom_params:
        if at not in counts:
            counts[at] = 0
        else:
            counts[at] += 1
        i = counts[at]

        # Get the corresponding atom block
        blocks = list()
        for _at, _params in atom_blocks:
            _at = "".join([i for i in _at if not i.isdigit()])
            if _at == at:
                blocks.append(_params)

        idx = 0 if len(blocks) == 1 else i
        at_params = blocks[idx].copy()
        at_params.update(params)
        at_params["symbol"] = at
        atoms.append(at_params)

    return atoms


class Atom:
    def __init__(
        self,
        symbol="",
        iq=1,
        it=1,
        ita=1,
        iz=0,
        conc=1.0,
        sms=1.0,
        sws=1.0,
        wswst=1.0,
        qtr=0.0,
        splt=0.0,
        fix="N",
        nz=None,
        norb=0,
        ion=0,
        config="",
        n=None,
        kappa=None,
        occup=None,
        valen=None,
    ):
        element = elements.get(symbol, {})
        iz = iz or element.get("iz", 0)
        nz = nz if nz is not None else iz
        ion = ion or element.get("ion", 0)
        norb = norb or element.get("norb", 0)
        config = config or element.get("config", "")
        n = n or element.get("n", [])
        kappa = kappa or element.get("kappa", [])
        occup = occup or element.get("occup", [])
        valen = valen or element.get("valen", [])

        self.symbol = symbol
        self.iz = iz
        self.norb = norb
        self.ion = ion
        self.config = config

        self.iq = iq  # Atomic site
        self.it = it  # Sublattice
        self.ita = ita  # Types of atoms occupying a given atomic site (in CPA)
        self.nz = nz  # Atomic number (same as IZ?)
        self.conc = conc  # Concentration of a type of atom in a given atomic site
        self.sms = sms  # Size of the local muffin-tin zero in units of S
        self.sws = sws  # Size of the potential spheres in units of WS
        self.wswst = wswst  # Size of atomic spheres in ASA in units of W.-S. spheres
        self.qtr = qtr  # Initial charge transfer
        self.splt = splt  # Initial magnetic moment
        self.fix = fix  # Fixed to the value of SPLT (Y) or it is not fixed (N) AFM=m

        self.u = list()
        self.j = list()

        self.n = np.array(n) if n else np.zeros(self.norb, dtype=np.int64)
        self.kappa = np.array(kappa) if n else np.zeros(self.norb, dtype=np.int64)
        self.occup = np.array(occup) if n else np.zeros(self.norb, dtype=np.int64)
        self.valen = np.array(valen) if n else np.zeros(self.norb, dtype=np.int64)

    def __repr__(self):
        cls = self.__class__.__name__
        s = (
            f"<{cls}({self.symbol} IQ={self.iq} It={self.it} Ita={self.ita}"
            f" Iz={self.iz} Norb={self.norb} Ion={self.ion} Config={self.config})>"
        )
        return s

    def __str__(self):
        cls = self.__class__.__name__
        s = (
            f"{cls}({self.symbol:<2} IQ={self.iq:<2} It={self.it:<2} Ita={self.ita:<2}"
            f" Iz={self.iz:<2} Norb={self.norb:<2} Ion={self.ion} Config={self.config})"
        )
        return s

    def set_n(self, a):
        self.n[:] = a

    def set_kappa(self, a):
        self.kappa[:] = a

    def set_occup(self, a):
        self.occup[:] = a

    def set_valen(self, a):
        self.valen[:] = a

    def set_num_valen(self, n):
        self.valen[:] = 0
        self.valen[-n:] = 1

    def init_dmft_energies(self, u: list = None, j: list = None):
        u = [0.0, 0.0, 0.0, 0.0] if u is None else u
        j = [0.0, 0.0, 0.0, 0.0] if j is None else j
        assert len(u) == len(j)
        self.u = list(u)
        self.j = list(j)

    def clear_dmft_energies(self):
        self.u.clear()
        self.j.clear()

    def get_u(self, item):
        if isinstance(item, str):
            item = ORBITALS.index(item)
        return self.u[item]

    def get_j(self, item):
        if isinstance(item, str):
            item = ORBITALS.index(item)
        return self.j[item]

    def set_u(self, item, u):
        if isinstance(item, str):
            item = ORBITALS.index(item)
        self.u[item] = u

    def set_j(self, item, j):
        if isinstance(item, str):
            item = ORBITALS.index(item)
        self.j[item] = j

    def update(self, *args, **kwargs):
        data = dict(*args, **kwargs)
        if "n" in data:
            self.n = np.array(data.pop("n"), dtype=np.int64)
        if "kappa" in data:
            self.kappa = np.array(data.pop("kappa"), dtype=np.int64)
        if "occup" in data:
            self.occup = np.array(data.pop("occup"), dtype=np.int64)
        if "valen" in data:
            self.valen = np.array(data.pop("valen"), dtype=np.int64)
        for k, v in data.items():
            self.__setattr__(k, v)

    def to_dict(self):
        data = self.__dict__.copy()
        return data

    @classmethod
    def from_dict(cls, data):
        self = cls()
        self.update(data)
        if self.norb:
            assert len(self.n) == self.norb
            assert len(self.kappa) == self.norb
            assert len(self.occup) == self.norb
            assert len(self.valen) == self.norb
            assert len(self.u) == len(self.j)
        else:
            assert len(self.n) == 1
            assert len(self.kappa) == 1
            assert len(self.occup) == 1
            assert len(self.valen) == 1
        return self

    def diff(self, other, exclude=None):
        return dict_diff(self.to_dict(), other.to_dict(), exclude)

    def __eq__(self, other):
        return self.to_dict() == other.to_dict()


def _frmt_atom_keys(key, iq, it, ita):
    keystrs = []
    if key is not None:
        keystrs.append(f"key={key}")
    if iq is not None:
        keystrs.append(f"iq={iq}")
    if it is not None:
        keystrs.append(f"it={it}")
    if ita is not None:
        keystrs.append(f"ita={ita}")
    return ", ".join(keystrs)


class KgrnFile(EmtoFile):
    """KGRN input file."""

    extension = ".dat"

    def __init__(self, path=None, dmft: bool = None, ga: bool = None, **kwargs):
        super().__init__(path)
        self._dmft = dmft
        self._ga = ga

        self.jobnam = "kgrn"
        self.header = ""  # first line after KGRN (usually date in the format %d %b %y)
        self.strt = "A"  # A: start from scratch, B: Resume, N: Reuse kmesh
        self.msgl = 1  # Level of printing
        self.expan = "S"  # Expansion mode: single (S), double (D) or modified (M)
        self.fcd = "Y"  # Y if full charge density is calculated, N if not
        self.func = "SCA"  # SCA (spherical cell approx), ASA (atomic sphere approx)

        # Input/Output directories
        self.for001 = ""
        self.for001_2 = ""
        self.for004 = ""
        self.dir002 = "pot/"
        self.dir003 = "pot/"
        self.dir006 = ""
        self.dir009 = "pot/"
        self.dir010 = "chd/"
        self.dir011 = ""

        # Line after path variables, eg 'Self-consistent KKR calculation for {jobnam}'
        self.comment = "Self-consistent KKR calculation for {jobnam}"

        self.niter = 50  # Max number of iterations in the main self-consistent DFT loop
        self.nlin = 30  # Max number of iterations in the Dyson equation.
        self.nprn = 0  # Determines what/how much will be printed in the output file
        self.ncpa = 7  # Maximum number of iterations in CPA loop.
        self.nt = 1  # Total number of sublattices: see atoms!
        self.mnta = 0  # Number of atomic species in an atomic site: see atoms!
        self.mode = "3D"  # 2D or 3D
        self.frc = "N"  # Whether (Y) or not (N) forces should be calculated.
        self.dos = "D"  # DOS (D), Fermi surface (F) or neither (N)
        self.ops = "N"  # Not used!
        self.afm = "P"  # Ferro (F), non (P), fixed total (M) or fixed spins (m)
        self.crt = "M"  # M/m for metallic systems, I for insulator
        self.lmaxh = 8  # Number of orbitals in the full charge density calculation.
        self.lmaxt = 4  # Number of tails.
        self.nfi = 31  # Number of mesh points in the Gaussian integration method.
        self.fixg = 2  # Determines the way the so muffin-tin zero V0 is calculated
        self.shf = 0  # Point where the stress tensors are calculated.
        self.sofc = "N"  # Y= soft core, N=frozen core, Z=all-electron frozen core
        self.kmsh = "G"  # k-mesh generation algorithm.
        self.ibz = 3  # Brillouin zone (same number as LAT in BMDL and KSTR)
        self.nkx, self.nky, self.nkz = 0, 0, 0  # The number of k points.
        self.fbz = "N"  # Y if the full Brillouin zone should be used, N if not.
        self.kmsh2 = "G"  # Not used!
        self.ibz2 = 1
        self.nkx2, self.nky2, self.nkz2 = 0, 0, 0
        self.zmsh = "C"  # Shape of the z-mesh along which the GF is calculated.
        self.nz1 = 16  # Ceiling for the number of points in z-mesh NZ=min(NZ+2,NZ1)
        self.nz2 = 16  # Parameter for the Fermi function integral
        self.nz3 = 8  # Parameter for the Fermi function integral
        self.nres = 4  # Parameter for the Fermi function integral
        self.nzd = 200  # Number of energy points in the DOS calculation
        self.depth = 1.0  # Width of z-mesh
        self.imagz = 0.01  # Imaginary broadening parameter
        self.eps = 0.2  # Parameter related to the energy integral.
        self.elim = -1.0  # Defines the crossing point of the z-mesh and the real axis
        # Mixing parameters for the DOS simple linear mixing scheme and Fermi energy
        self.amix, self.efmix = 0.01, 1.0
        self.vmtz = 0.0  # Shift of the Muffin-tin zero energy (FIXG=3)
        self.mmom = 0.0  # Total magnetic moment for fixed spin calculation (AFM=M)
        # Convergence criterions of the energy, Fermi energy and CPA loop
        self.tole, self.tolef, self.tolcpa = 1e-7, 1e-7, 1e-6
        self.tfermi = 500  # Fermi temperature in Kelvins. Used when ZMSH=F
        self.sws = 0.0  # Average Wigner-Seitz radius, volume of the unit cell (bohr)
        self.nsws = 1  # *Not used!*
        self.dsws = 0.05  # *Not used!*
        self.alpcpa = 0.602  # The screened impurity model parameter of CPA
        self.efgs = 0.0  # Initial guess for Fermi level
        self.hx = 0.1  # The distance between points on the linear path
        self.nx = 5  # The number of trial points on the z-mesh linear path
        self.nz0 = 16  # Initial number of z-mesh points.
        self.stmp = "N"  # Y: tmp disk storage (DIR011), A: RAM, N: No tmp storage

        self.iex = 4  # Determines which exchange-correlation functional to use (4=LDA)
        self.np = 251  # number of radial grid points in the Poisson equation solver.
        # Number of times the atomic orbital energies are adjusted in the Dirac solver
        self.nes = 15
        self.dirac_niter = 100  # The maximum number of iterations in the Dirac solver
        # If IWAT = 1, the potential DV(I) is corrected either with a term ION/RWAT or
        # ION/DR(I), where ION is the number of electrons of a given atomic species.
        # DR(I) is the point I of the radial mesh. If IWAT = 0,
        # this correction is not performed.
        self.iwat = 0
        self.nprna = 0  # Determines how much information will be printed in the file
        self.vmix = 0.3  # Mixing parameter of the potential.
        self.rwat = 3.5  # The radius of the Watson sphere.
        # Controls the distribution of points in the atomic radial mesh
        self.rmax = 20.0
        # A step size controlling the number of points in the atomic radial mesh.
        # Small values of DX lead to high number of points and vice versa.
        self.dx = 0.03
        self.dr1 = 0.002
        # Convergence criteria in the Poisson equation and the orbital Dirac equations
        self.test = self.teste = self.testy = self.testv = 1e-12

        self.avg = "A"  # Average type, A: arithmetic (CPA), G: geometric (TMT)

        # Atoms
        self.atoms = list()

        # Optional DMFT parameters
        self.for007 = ""
        self.nom = 1024
        self.nomi = 40
        self.dc = 1
        self.ttt = 400.0
        self.smix = 0.5
        self.solver = "uppsalasolver"

        self.load(missing_ok=True)
        if kwargs:
            self.update(kwargs)

    @property
    def is_dmft(self) -> bool:
        return self._dmft

    def init(self, jobname: str, header: str = "", **kwargs) -> None:
        self.jobnam = jobname
        self.set_header(header)
        self.update(kwargs)

    def force_dmft(self, dmft: bool = True) -> None:
        self._dmft = dmft

    def force_ga(self, ga: bool = True) -> None:
        if ga and not self._dmft:
            raise ValueError("DMFT must be enabled to use GA")
        self._ga = ga

    def aux_dirs(self) -> List[str]:
        d = self.dir002, self.dir003, self.dir006, self.dir009, self.dir010, self.dir011
        return [path for path in d if path]

    def set_header(self, header: str = "", date_frmt: str = "%d %b %y") -> None:
        self.header = (header + " " + datetime.now().strftime(date_frmt)).strip()

    def get_atom_symbols(self, include_empty: bool = True) -> List[str]:
        atoms = list()
        for at in self.atoms:
            sym = at.symbol
            if include_empty or sym not in ("E", "Va"):
                if sym not in atoms:
                    atoms.append(sym)
        return atoms

    def get_atoms(
        self,
        key: Union[int, str] = None,
        iq: int = None,
        it: int = None,
        ita: int = None,
    ) -> List[Atom]:
        # Integer key
        if isinstance(key, int):
            if any([iq, it, ita]):
                raise ValueError(
                    "Integer key and iq, it or ita cannot be used together"
                )
            try:
                return [self.atoms[key]]
            except IndexError:
                raise KeyError(f"Atom with index {key} not found")
        # Other cases
        atoms = list()
        for at in self.atoms:
            if key is not None and at.symbol.lower() != key.lower():
                continue
            if iq is not None and at.iq != iq:
                continue
            if it is not None and at.it != it:
                continue
            if ita is not None and at.ita != ita:
                continue
            atoms.append(at)
        if not atoms:
            keystr = _frmt_atom_keys(key, iq, it, ita)
            raise KeyError(f"No atoms found with {keystr}")
        return atoms

    def get_atom(
        self,
        key: Union[int, str] = None,
        iq: int = None,
        it: int = None,
        ita: int = None,
    ) -> Atom:
        atoms = self.get_atoms(key, iq, it, ita)
        if len(atoms) > 1:
            keystr = _frmt_atom_keys(key, iq, it, ita)
            raise KeyError(f"Multiple atoms found with {keystr}: {atoms}")
        return atoms[0]

    def add_atom(self, atom: Union[str, Atom], **kwargs) -> Atom:
        if isinstance(atom, str):
            atom = Atom(atom)
        if kwargs:
            atom.update(**kwargs)
        self.atoms.append(atom)
        return atom

    def set_concentrations(self, cc: List[float]) -> None:
        assert len(cc) == len(self.atoms)
        for atom, conc in zip(self.atoms, cc):
            atom.conc = conc

    def get_concentrations(self) -> List[float]:
        return [atom.conc for atom in self.atoms]

    def get_concentration(self, key: Union[int, str]) -> float:
        try:
            atoms = self.get_atoms(key)
            concs = np.unique([a.conc for a in atoms])
            if len(concs) > 1:
                raise ValueError(f"Multiple concentrations found: {concs}")
            return float(concs[0])
        except KeyError:
            return 0.0

    def set_uj(self, atom, u, j, idx=None):
        atoms = list()
        for at in self.get_atoms(atom):
            if idx is None:
                at.u = u
                at.j = j
            else:
                at.u[idx] = u
                at.j[idx] = j
            atoms.append(at)
        return atoms

    def get_max_it(self) -> int:
        return max(at.it for at in self.atoms)

    def get_max_ita(self, it: int = None) -> int:
        if it is None:
            atoms = self.atoms
        else:
            atoms = self.get_atoms(it=it)
        return max(at.ita for at in atoms)

    def update_mnta(self) -> None:
        self.mnta = max(atom.ita for atom in self.atoms)

    def init_dmft_energies(self, u: List[float] = None, j: List[float] = None) -> None:
        for atom in self.atoms:
            atom.init_dmft_energies(u=u, j=j)

    def clear_dmft_energies(self) -> None:
        for atom in self.atoms:
            atom.clear_dmft_energies()

    def set_kstr_path(self, path: Union[str, Path]) -> None:
        path = str(path)
        assert path.endswith(".tfh")
        self.for001 = path

    def set_kstr2_path(self, path: Union[str, Path]) -> None:
        path = str(path)
        assert path.endswith(".tfh")
        self.for001_2 = path

    def set_bmdl_path(self, path: Union[str, Path]) -> None:
        path = str(path)
        assert path.endswith(".mdl")
        self.for004 = path

    def set_pot_path(self, path: Union[str, Path]) -> None:
        path = str(path)
        self.dir002 = path
        self.dir003 = path
        self.dir009 = path

    def set_chd_path(self, path: Union[str, Path]) -> None:
        path = str(path)
        self.dir010 = path

    def set_tmp_path(self, path: Union[str, Path]) -> None:
        path = str(path)
        self.dir011 = path

    def update_paths(
        self,
        kstr: Union[str, Path],
        bmdl: Union[str, Path],
        kstr2: Union[str, Path] = "",
        pot: Union[str, Path] = "pot/",
        chd: Union[str, Path] = "chd/",
        tmp: Union[str, Path] = "",
    ) -> None:
        kstr = str(kstr)
        bmdl = str(bmdl)
        kstr2 = str(kstr2)
        pot = str(pot)
        chd = str(chd)
        tmp = str(tmp)
        assert kstr.endswith(".tfh")
        assert bmdl.endswith(".mdl")

        self.for001 = kstr
        if kstr2:
            assert kstr2.endswith(".tfh")
            self.for001_2 = kstr2
        else:
            self.for001_2 = ""
        self.for004 = bmdl
        self.dir002 = pot
        self.dir003 = pot
        self.dir006 = ""
        self.dir009 = pot
        self.dir010 = chd
        self.dir011 = tmp

    def update_from_config(
        self,
        kstr: Union[str, Path],
        bmdl: Union[str, Path],
        kstr2: Union[str, Path] = "",
        pot: Union[str, Path] = "pot/",
        chd: Union[str, Path] = "chd/",
        tmp: Union[str, Path] = "",
        conf: dict = None,
    ):
        kstr = str(kstr)
        bmdl = str(bmdl)
        kstr2 = str(kstr2)
        pot = str(pot)
        chd = str(chd)
        tmp = str(tmp)
        update_emto_paths(
            self, kstr, bmdl, kstr2=kstr2, pot=pot, chd=chd, tmp=tmp, conf=conf
        )

    def param_diff(self, other, exclude=None):
        if isinstance(other, KgrnFile):
            other = other.to_dict()
        elif not isinstance(other, dict):
            raise TypeError(f"other has to be a KgrnFile or dict, not {type(other)}")
        return dict_diff(self.to_dict(), other, exclude=exclude)

    def update(self, *args, **kwargs):
        data = dict(*args, **kwargs)
        for k, v in data.items():
            if not hasattr(self, k):
                raise KeyError(f"{k} is not a valid field of {self.__class__.__name__}")
            if isinstance(v, str):
                v = v.strip()
            self.__setattr__(k, v)

    def to_dict(self) -> Dict[str, Union[str, int, float]]:
        data = {k: v for k, v in self.__dict__.items() if not k.startswith("_")}
        data.pop("atoms", None)
        data.pop("path", None)
        if not self.is_dmft:
            dmft_keys = ["for007", "nom", "nomi", "dc", "ttt", "smix", "solver"]
            for key in dmft_keys:
                data.pop(key, None)
        return data

    def __getitem__(self, key: str) -> Union[str, int, float]:
        key = key.strip().lower()
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        return self.__getattribute__(key)

    def __setitem__(self, key: str, value: Union[str, int, float]) -> None:
        key = key.strip().lower()
        # Try to convert type
        dtype = type(self.__getattribute__(key))
        try:
            value = dtype(value)
        except ValueError:
            pass
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        self.__setattr__(key, value)

    def check(self) -> None:
        """Check if the input is consistent."""
        if self.jobnam is None:
            raise KGRNError("'jobnam' has to be given!")
        if self.strt not in ("A", "B", "N"):
            raise KGRNError(f"'strt' has to be 'A', 'B' or 'N', not '{self.strt}'!")
        if self.expan not in ("S", "D", "M"):
            raise KGRNError(f"'expan' has to be 'S', 'D' or 'M', not '{self.expan}'!")
        if self.expan == "M" and not self.for001_2:
            raise KGRNError("for001_2 has to be given when expan='M'!")
        if self.fcd not in ("Y", "N"):
            raise KGRNError(f"'fcd' has to be 'Y' or 'N', not '{self.fcd}'!")
        if self.func not in ("SCA", "ASA"):
            raise KGRNError(f"'func' has to be 'SCA' or 'ASA', not {self.func}!")
        if self.avg not in ("A", "G"):
            raise KGRNError(f"'avg' has to be 'A' or 'G', not {self.avg}!")
        if len(self.atoms) == 0:
            raise KGRNError("No atoms given!")
        # Check if atoms are consistent with parameters
        nt = max(atom.it for atom in self.atoms)
        if self.nt != nt:
            raise KGRNError(f"maximal iq={nt} does not match nt={self.nt}!")
        mnta = max(atom.ita for atom in self.atoms)
        if self.mnta != mnta:
            raise KGRNError(f"maximal ita={mnta} does not match mnta={self.mnta}!")

    # ----------------------------------------------------------------------------------

    def loads(self, text: str) -> "KgrnFile":
        # get template by file key (first 4 characters, KGRN or DMFT)
        fkey = text[:4]
        if fkey == "KGRN":
            tmplt_str = TEMPLATE_KGRN
            if self._dmft is None:
                self._dmft = False
        elif fkey == "DMFT":
            tmplt_str = TEMPLATE_DMFT
            if self._dmft is None:
                self._dmft = True
        else:
            raise KGRNReadError(f"Unknown file format: Invalid key '{fkey}'")

        if self._dmft and "AVG..=" in text:
            self._ga = True
            tmplt_str = TEMPLATE_DMFT_v2

        template = Template(tmplt_str, ignore_case=True)

        # Parse file contents
        try:
            data = template.parse(text.replace(".d", ".e"))
        except Exception as e:
            raise KGRNReadError(f"Failed to parse file: {self.path}") from e

        params = dict(data.copy())
        atomstr = params.pop("atoms")
        confstr = params.pop("atomconf")
        try:
            atoms = parse_atoms(atomstr, confstr)
        except Exception as e:
            raise KGRNReadError(f"Failed to parse atoms: {self.path}\n{atomstr}") from e

        self.update(params)
        self.atoms = [Atom.from_dict(at) for at in atoms]
        # self.check()
        return self

    def dumps(self) -> str:
        dmft = self._dmft

        # get template by file key (first 4 characters, KGRN or DMFT)
        if dmft:
            if self._ga:
                tmplt_str = TEMPLATE_DMFT_v2
            else:
                tmplt_str = TEMPLATE_DMFT
        else:
            tmplt_str = TEMPLATE_KGRN

        template = Template(tmplt_str, ignore_case=True)

        # Check if input is consistent and convert to dict
        try:
            self.check()
        except KGRNError as e:
            raise KGRNWriteError(f"Invalid input: {e}") from e
        params = self.to_dict()

        # Pre-format header and comment fields to allow for nested formatting
        params["header"] = self.header.format(**params)
        params["comment"] = self.comment.format(**params)

        # Format atom info
        atomstr = "Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT"
        if not self.atoms or self.atoms[0].fix:
            atomstr += " Fix"
        atomstr += "\n"
        atoms = [at.to_dict() for at in self.atoms]
        atomstr += "\n".join(format_atom_line(atom) for atom in atoms)
        atomconf = "\n".join(format_atom_block(atom) for atom in atoms)
        data = dict(params.copy())
        data["atoms"] = atomstr
        data["atomconf"] = atomconf
        try:
            return template.format(data).replace(".0e", ".d")
        except Exception as e:
            raise KGRNWriteError(f"Failed to format file: {self.path}") from e
