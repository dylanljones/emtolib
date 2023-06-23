# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

import re
from datetime import datetime
import logging
import numpy as np
from typing import Union
from ..common import EmtoFile, parse_params, elements
from ..ftmplt import Template

logger = logging.getLogger(__name__)

RE_BAND_SECTION = re.compile(r"Band: (.*?) lines")

TEMPLATE = """\
KGRN {header:>55}
JOBNAM={jobnam}
STRT..=  {strt} MSGL.=  {msgl:d} EXPAN.= {expan} FCD..=  {fcd} FUNC..= {func}
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
NITER.={niter:3d} NLIN.={nlin:3d} NPRN.=  {nprn:d} NCPA.={ncpa:3d} NT...={nt:3d} MNTA.={mnta:3d}
MODE..= {mode:2} FRC..=  {frc} DOS..=  {dos} OPS..=  {ops} AFM..=  {afm} CRT..=  {crt}
Lmaxh.={lmaxh:3d} Lmaxt={lmaxt:3d} NFI..={nfi:3d} FIXG.={fixg:3d} SHF..=  {shf:1d} SOFC.=  {sofc}
KMSH...= {kmsh} IBZ..={ibz:3d} NKX..={nkx:3d} NKY..={nky:3d} NKZ..={nkz:3d} FBZ..=  {fbz}
KMSH2..= {kmsh2} IBZ2.={ibz2:3d} NKX2.={nkx2:3d} NKY2.={nky2:3d} NKZ2.={nkz2:3d}
ZMSH...= {zmsh} NZ1..={nz1:3d} NZ2..={nz2:3d} NZ3..={nz3:3d} NRES.={nres:3d} NZD.={nzd:4d}
DEPTH..= {depth:6.3f} IMAGZ.={imagz:7.5f} EPS...={eps:7.5f} ELIM..= {elim:6.3f}
AMIX...= {amix:6.3f} EFMIX.= {efmix:6.3f} VMTZ..={vmtz:7.3f} MMOM..={mmom:7.3f}
TOLE...= {tole:7.1e} TOLEF.= {tolef:7.1e} TOLCPA= {tolcpa:7.1e} TFERMI= {tfermi:6.1f} (K)
SWS......={sws:10.7f} NSWS.={nsws:3d} DSWS..=   {dsws:4.2f} ALPCPA= {alpcpa:6.4f}
Setup: 2 + NQ*NS lines
EFGS...= {efgs:6.3f} HX....= {hx:6.3f} NX...= {nx:2d} NZ0..= {nz0:2d} STMP..= {stmp}
{atoms}
Atom:  4 lines + NT*NTA*6 lines
IEX...= {iex:2d} NP..={np:4d} NES..={nes:3d} NITER={dirac_niter:3d} IWAT.={iwat:3d} NPRNA={nprna:3d}
VMIX.....=  {vmix:8.6f} RWAT....=  {rwat:8.6f} RMAX....={rmax:10.6f}
DX.......=  {dx:8.6f} DR1.....=  {dr1:8.6f} TEST....=  {test:8.2E}
TESTE....=  {teste:8.2E} TESTY...=  {testy:8.2E} TESTV...=  {testv:8.2E}
{atomconf}
"""  # noqa

ATOM_COLUMNS = tuple("iq it ita nz conc sms sws wswst qtr splt fix".split())
ATCONF_TEMPLATE = "Iz= {iz:3d} Norb={norb:3d} Ion=  {ion:d} Config= {config}"
ATLINE_TEMPLATE = (
    "{symbol:2}    {iq:3d} {it:2d} {ita:2d}  {nz:2d}  {conc:5.3f}  "
    "{sms:5.3f}  {sws:5.3f}  {wswst:5.3f} {qtr:4.1f} {splt:4.1f}   {fix}"
)

ORBITALS = ["s", "p", "d", "f", "g", "h"]


class KGRNError(Exception):
    pass


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
        raise KGRNError(f"Invalid atom header: {header}")
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
        n = n or [int(x) for x in element.get("n", "").split() if x.strip()]
        kappa = kappa or [int(x) for x in element.get("kappa", "").split() if x.strip()]
        occup = occup or [int(x) for x in element.get("occup", "").split() if x.strip()]
        valen = valen or [int(x) for x in element.get("valen", "").split() if x.strip()]

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
            f"<{cls}({self.symbol} It={self.it} Ita={self.ita} Iz={self.iz} "
            f"Norb={self.norb} Ion={self.ion} Config={self.config})>"
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


class EmtoKgrnFile(EmtoFile):
    """KGRN input file."""

    extension = ".dat"
    template = Template(TEMPLATE, ignore_case=True)

    def __init__(self, path=None, **kwargs):
        super().__init__(path)

        self.jobnam = "kgrn"
        self.header = ""  # first line after KGRN (usually date in the format %d %b %y)
        self.strt = "A"  # A: start from scratch, B: Resume, N: Reuse kmesh
        self.msgl = 1  # Level of printing
        self.expan = "S"  # Expansion mode: single (S), double (D) or modified (M)
        self.fcd = "Y"  # Y if full charge density is calculated, N if not
        self.func = "SCA"  # SCA (spherical cell approx), ASA (atomic sphere approx)

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
        self.alpcpa = 0.6  # The screened impurity model parameter of CPA
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

        # Input/Output directories
        self.for001 = ""
        self.for001_2 = None
        self.for004 = ""
        self.dir002 = "pot/"
        self.dir003 = "pot/"
        self.dir006 = ""
        self.dir009 = "pot/"
        self.dir010 = "chd/"
        self.dir011 = "/tmp/"

        self.atoms = list()

        if self.path.is_file():
            self.load()

        if kwargs:
            self.update(kwargs)

    def aux_dirs(self):
        d = self.dir002, self.dir003, self.dir006, self.dir009, self.dir010, self.dir011
        return [path for path in d if path]

    def set_header(self, header="", date_frmt="%d %b %y"):
        self.header = (header + " " + datetime.now().strftime(date_frmt)).strip()

    def get_atom_symbols(self, include_empty=True):
        atoms = set(
            at.symbol
            for at in self.atoms
            if include_empty or at.symbol not in ("E", "Va")
        )
        return atoms

    def get_atom(self, key) -> Atom:
        if isinstance(key, int):
            return self.atoms[key]
        # Check if index in key
        if not key.isalpha():
            # split idx from jey
            idx = re.findall(r"\d+", key)[0]
            key = key.replace(idx, "").strip()
            idx = int(idx.strip())
        else:
            idx = 0
        res = list()
        for at in self.atoms:
            if at.symbol == key:
                res.append(at)
        return res[idx]

    def add_atom(self, atom: Union[str, Atom]):
        if isinstance(atom, str):
            atom = Atom(atom)
        self.atoms.append(atom)
        return atom

    def set_concentrations(self, cc):
        assert len(cc) == len(self.atoms)
        for atom, conc in zip(self.atoms, cc):
            atom.conc = conc

    def get_concentrations(self):
        return [atom.conc for atom in self.atoms]

    def get_concentration(self, key):
        try:
            return self.get_atom(key).conc
        except IndexError:
            return 0.0

    def param_diff(self, other, exclude=None):
        if exclude is not None and isinstance(exclude, str):
            exclude = (exclude,)
        d1 = self.to_dict()
        d2 = other.to_dict()
        diffset = set(d1.items()) ^ set(d2.items())
        diff = dict()
        for key in dict(diffset).keys():
            if exclude is None or key not in exclude:
                diff[key] = (d1[key], d2[key])
        return diff

    def update(self, *args, **kwargs):
        data = dict(*args, **kwargs)
        for k, v in data.items():
            if not hasattr(self, k):
                raise KeyError(f"{k} is not a valid field of {self.__class__.__name__}")
            if isinstance(v, str):
                v = v.strip()
            self.__setattr__(k, v)

    def check(self):
        """Check if the input is consistent."""
        if self.jobnam is None:
            raise KGRNError("'jobnam' has to be given!")
        if self.strt not in ("A", "B", "N"):
            raise KGRNError("'strt' has to be 'A', 'B' or 'N'!")
        if self.expan not in ("S", "D", "M"):
            raise KGRNError("'expan' has to be 'S', 'D' or 'M'!")
        if self.fcd not in ("Y", "N"):
            raise KGRNError("'fcd' has to be 'Y' or 'N'!")
        if self.func not in ("SCA", "ASA"):
            raise KGRNError("'func' has to be 'SCA' or 'ASA'!")

    # ----------------------------------------------------------------------------------

    def to_dict(self):
        data = {k: v for k, v in self.__dict__.items() if not k.startswith("_")}
        data.pop("atoms", None)
        data.pop("path", None)
        if data["for001_2"] is None:
            data["for001_2"] = data["for001"]
        return data

    def __getitem__(self, key):
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        return self.__getattribute__(key)

    def __setitem__(self, key, value):
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        self.__setattr__(key, value)

    def loads(self, text):
        data = self.template.parse(text.replace(".d", ".e"))
        params = dict(data.copy())
        atomstr = params.pop("atoms")
        confstr = params.pop("atomconf")
        try:
            atoms = parse_atoms(atomstr, confstr)
        except Exception as e:
            raise KGRNError(f"Failed to parse atoms: {self.path}\n{atomstr}") from e

        self.update(params)
        self.atoms = [Atom.from_dict(at) for at in atoms]
        self.check()
        return self

    def dumps(self):
        self.check()
        params = self.to_dict()

        # Pre-format header and comment fields to allow for nested formatting
        params["header"] = self.header.format(**params)
        params["comment"] = self.comment.format(**params)

        # format atom info
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

        return self.template.format(data).replace(".0e", ".d")
