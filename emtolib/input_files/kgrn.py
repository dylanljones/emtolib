# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

import re
from datetime import datetime
import logging
import numpy as np
from ..common import EmtoFile, parse_params, elements

logger = logging.getLogger(__name__)

RE_BAND_SECTION = re.compile(r"Band: (.*?) lines")

TEMPLATE = """\
KGRN                                               {date:%d %b %y}
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
Self-consistent KKR calculation for {jobnam}
Band: 10 lines
NITER.={niter:3d} NLIN.={nlin:3d} NPRN.=  {nprn:d} NCPA.={ncpa:3d} NT...= {nt:2d} MNTA.= {mnta:2d}
MODE..= {mode:2} FRC..=  {frc} DOS..=  {dos} OPS..=  {ops} AFM..=  {afm} CRT..=  {crt}
Lmaxh.= {lmaxh:2d} Lmaxt= {lmaxt:2d} NFI..={nfi:3d} FIXG.= {fixg:2d} SHF..=  {shf:1d} SOFC.=  {sofc}
KMSH...= {kmsh} IBZ..= {ibz:2d} NKX..= {nkx:2d} NKY..= {nky:2d} NKZ..= {nkz:2d} FBZ..=  {fbz}
KMSH2..= {kmsh2} IBZ2.={ibz2:3d} NKX2.={nkx2:3d} NKY2.={nky2:3d} NKZ2.={nkz2:3d}
ZMSH...= {zmsh} NZ1..= {nz1:2d} NZ2..={nz2:3d} NZ3..={nz3:3d} NRES.={nres:3d} NZD.={nzd:4d}
DEPTH..= {depth:6.3f} IMAGZ.={imagz:7.5f} EPS...={eps:7.5f} ELIM..= {elim:6.3f}
AMIX...= {amix:6.3f} EFMIX.= {efmix:6.3f} VMTZ..={vmtz:7.3f} MMOM..={mmom:7.3f}
TOLE...= {tole:7.1e} TOLEF.= {tolef:7.1e} TOLCPA= {tolcpa:7.1e} TFERMI= {tfermi:6.1f} (K)
SWS......={sws:10.7f} NSWS.={nsws:3d} DSWS..=   {dsws:4.2f} ALPCPA= {alpcpa:6.4f}
Setup: 2 + NQ*NS lines
EFGS...= {efgs:6.3f} HX....= {hx:6.3f} NX...= {nx:2d} NZ0..= {nz0:2d} STMP..= {stmp}
Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT Fix
{atoms}
Atom:  4 lines + NT*NTA*6 lines
IEX...= {iex:2d} NP..={np:4d} NES..={nes:3d} NITER={dirac_niter:3d} IWAT.={iwat:3d} NPRNA={nprna:3d}
VMIX.....=  {vmix:8.6f} RWAT....=  {rwat:8.6f} RMAX....={rmax:10.6f}
DX.......=  {dx:8.6f} DR1.....=  {dr1:8.6f} TEST....=  {test:8.2E}
TESTE....=  {teste:8.2E} TESTY...=  {testy:8.2E} TESTV...=  {testv:8.2E}
{atomconf}
"""  # noqa

ATOM_COLUMNS = "iq it ita nz conc sms sws wswst qtr splt fix".split()
ATCONF_TEMPLATE = "Iz= {iz:3d} Norb={norb:3d} Ion=  {ion:d} Config= {config}"
ATLINE_TEMPLATE = (
    "{symbol:2}    {iq:3d} {it:2d} {ita:2d}  {nz:2d}  {conc:5.3f}  "
    "{sms:5.3f}  {sws:5.3f}  {wswst:5.3f} {qtr:4.1f}{splt:5.2f}   {fix}"
)

ORBITALS = ["s", "p", "d", "f", "g", "h"]


def format_atom_line(params):
    line = ATLINE_TEMPLATE.format(**params)
    uj = params["u"] + params["j"]
    if uj:
        line += " " + " ".join(f"{x:4.2f}" for x in uj)
    return line


def format_atom_block(params):
    lines = list()
    lines.append(params["symbol"])
    lines.append(ATCONF_TEMPLATE.format(**params))
    lines.append("n     " + " ".join(f"{x:2d}" for x in params["n"]))
    lines.append("Kappa " + " ".join(f"{x:2d}" for x in params["kappa"]))
    lines.append("Occup " + " ".join(f"{x:2d}" for x in params["occup"]))
    lines.append("Valen " + " ".join(f"{x:2d}" for x in params["valen"]))
    return "\n".join(lines)


def parse_atom_line(line):
    columns = ATOM_COLUMNS
    ncols = len(columns)
    parts = line.split()
    at = parts.pop(0)
    values = parts[:ncols]
    for i in range(4):
        values[i] = int(values[i])
    for i in range(4, ncols - 1):
        values[i] = float(values[i])
    params = dict(zip(columns, values))
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
    atom_params = [parse_atom_line(line) for line in atomstr.splitlines()]

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
        for _at, at_params in atom_blocks:
            _at = "".join([i for i in _at if not i.isdigit()])
            if _at == at:
                blocks.append(at_params)

        idx = 0 if len(blocks) == 1 else i
        at_params = blocks[idx]
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
        nz=None,
        conc=1.0,
        sms=1.0,
        sws=1.0,
        wswst=1.0,
        qtr=0.0,
        splt=0.0,
        fix="N",
        iz=None,
        norb=0,
        ion=0,
        config=None,
    ):
        element = elements.get(symbol, {})
        config = config or element.get("econf", "")
        nz = nz or element.get("iz", 0)
        iz = iz or nz

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

        self.n = np.zeros(self.norb, dtype=np.int64)
        self.kappa = np.zeros(self.norb, dtype=np.int64)
        self.occup = np.zeros(self.norb, dtype=np.int64)
        self.valen = np.zeros(self.norb, dtype=np.int64)

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
        assert len(self.n) == self.norb
        assert len(self.kappa) == self.norb
        assert len(self.occup) == self.norb
        assert len(self.valen) == self.norb
        assert len(self.u) == len(self.j)
        return self


class EmtoKgrnFile(EmtoFile):

    extension = ".dat"

    def __init__(self, path=None, update_date=True, **kwargs):
        super().__init__(path)
        self._update_date = update_date

        self.jobnam = "kgrn"
        self.date = datetime.now()
        self.strt = "A"  # A: start from scratch, B: Resume, N: Reuse kmesh
        self.msgl = 1  # Level of printing
        self.expan = "S"  # Expansion mode: single (S), double (D) or modified (M)
        self.fcd = "Y"  # Y if full charge density is calculated, N if not
        self.func = "SCA"  # SCA (spherical cell approx), ASA (atomic sphere approx)

        self.niter = 50
        self.nlin = 30
        self.nprn = 0
        self.ncpa = 7
        self.nt = 1  # Total number of sublattices: see atoms!
        self.mnta = 2  # Number of atomic species in an atomic site: see atoms!
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

        self.iex = 4
        self.np = 251
        self.nes = 15
        self.dirac_niter = 100
        self.iwat = 0
        self.nprna = 0
        self.vmix = 0.3
        self.rwat = 3.5
        self.rmax = 20.0
        self.dx = 0.03
        self.dr1 = 0.002
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

        if self.path.exists():
            self.load()

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

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)

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

    def param_diff(self, other):
        d1 = self.to_dict()
        d2 = other.to_dict()
        diffset = set(d1.items()) ^ set(d2.items())
        return {key: (d1[key], d2[key]) for key in dict(diffset).keys()}

    # ----------------------------------------------------------------------------------

    def to_dict(self):
        data = {k: v for k, v in self.__dict__.items() if not k.startswith("_")}
        data.pop("atoms", None)
        data.pop("path", None)
        return data

    def __getitem__(self, key):
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        return self.__getattribute__(key)

    def __setitem__(self, key, value):
        if not hasattr(self, key):
            raise KeyError(f"{key} is not a valid field of {self.__class__.__name__}")
        self.__setattr__(key, value)

    def update(self, *args, **kwargs):
        data = dict(*args, **kwargs)
        for k, v in data.items():
            if not hasattr(self, k):
                raise KeyError(f"{k} is not a valid field of {self.__class__.__name__}")
            self.__setattr__(k, v)

    def loads(self, data: str) -> None:
        params = dict()
        lines = data.splitlines(keepends=False)
        prog, date = lines.pop(0).split(" ", maxsplit=1)
        assert prog == "KGRN"
        self.date = datetime.strptime(date.strip(), "%d %b %y")
        self.jobnam = lines.pop(0).split("=")[1].strip()
        params.update(parse_params(lines.pop(0)))
        self.for001 = lines.pop(0).split("=")[1].strip()
        self.for001_2 = lines.pop(0).split("=")[1].strip()
        self.dir002 = lines.pop(0).split("=")[1].strip()
        self.dir003 = lines.pop(0).split("=")[1].strip()
        self.for004 = lines.pop(0).split("=")[1].strip()
        self.dir006 = lines.pop(0).split("=")[1].strip()
        self.dir009 = lines.pop(0).split("=")[1].strip()
        self.dir010 = lines.pop(0).split("=")[1].strip()
        self.dir011 = lines.pop(0).split("=")[1].strip()
        assert lines.pop(0).startswith("Self-consistent KKR calculation for")
        # Band section
        nlines = int(RE_BAND_SECTION.match(lines.pop(0)).group(1))
        params.update(parse_params(" ".join(lines[:nlines]).replace(" (K)", "")))
        # Setup section
        lines = lines[nlines:]
        line = lines.pop(0)
        assert line.startswith("Setup: ")
        nlines = int(re.findall(r"\d", line)[0]) - 1
        params.update(parse_params(" ".join(lines[:nlines])))

        # Atom lines
        lines = lines[nlines:]
        columns = lines.pop(0).split()[1:12]  # header line
        columns = [c.strip().lower().replace("(", "").replace(")", "") for c in columns]
        assert columns == ATOM_COLUMNS
        atomstr_lines = list()
        while not lines[0].startswith("Atom:"):
            atomstr_lines.append(lines.pop(0))
        atomstr = "\n".join(atomstr_lines)

        # Atom section
        line = lines.pop(0)
        assert line.startswith("Atom: ")
        nlines, nblock = (int(x) for x in re.findall(r"\d", line))
        atparams = parse_params(" ".join(lines[:nlines]))
        atparams["ATNITER"] = atparams["NITER"]
        del atparams["NITER"]
        params.update(atparams)

        # Atom blocks
        atomconfstr = "\n".join(lines[nlines:])

        params = {k.upper(): v for k, v in params.items()}
        self.strt = params["STRT"]
        self.msgl = int(params["MSGL"])
        self.expan = params["EXPAN"]
        self.fcd = params["FCD"]
        self.func = params["FUNC"]
        self.niter = int(params["NITER"])
        self.nlin = int(params["NLIN"])
        self.nprn = int(params["NPRN"])
        self.ncpa = int(params["NCPA"])
        self.nt = int(params["NT"])
        self.mnta = int(params["MNTA"])
        self.mode = params["MODE"]
        self.frc = params["FRC"]
        self.dos = params["DOS"]
        self.ops = params["OPS"]
        self.afm = params["AFM"]
        self.crt = params["CRT"]
        self.lmaxh = int(params["LMAXH"])
        self.lmaxt = int(params["LMAXT"])
        self.nfi = int(params["NFI"])
        self.fixg = int(params["FIXG"])
        self.shf = int(params["SHF"])
        self.sofc = params["SOFC"]
        self.kmsh = params["KMSH"]
        self.ibz = int(params["IBZ"])
        self.nkx, self.nky, self.nkz = (int(params[k]) for k in ["NKX", "NKY", "NKZ"])
        self.fbz = params["FBZ"]
        self.kmsh2 = params.get("KMSH2", "G")
        self.nkx2 = int(params.get("NKX2", "0"))
        self.nky2 = int(params.get("NKY2", "0"))
        self.nkz2 = int(params.get("NKZ2", "0"))
        self.zmsh = params["ZMSH"]
        self.nz1, self.nz2, self.nz3 = (int(params[k]) for k in ["NZ1", "NZ2", "NZ3"])
        self.nres = int(params["NRES"])
        self.nzd = int(params["NZD"])
        self.depth = float(params["DEPTH"])
        self.imagz = float(params["IMAGZ"])
        self.eps = float(params["EPS"])
        self.elim = float(params["ELIM"])
        self.amix = float(params["AMIX"])
        self.efmix = float(params["EFMIX"])
        self.vmtz = float(params["VMTZ"])
        self.mmom = float(params["MMOM"])
        self.tole = float(params["TOLE"].replace("d", "e"))
        self.tolef = float(params["TOLEF"].replace("d", "e"))
        self.tolcpa = float(params["TOLCPA"].replace("d", "e"))
        self.tfermi = float(params["TFERMI"])
        self.sws = float(params["SWS"])
        self.nsws = int(params["NSWS"])
        self.dsws = float(params["DSWS"])
        self.alpcpa = float(params["ALPCPA"])
        self.efgs = float(params["EFGS"])
        self.hx = float(params["HX"])
        self.nx = int(params["NX"])
        self.nz0 = int(params["NZ0"])
        self.stmp = params["STMP"]
        self.iex = int(params["IEX"])
        self.np = int(params["NP"])
        self.nes = int(params["NES"])
        self.dirac_niter = int(params["ATNITER"])
        self.iwat = int(params["IWAT"])
        self.nprna = int(params["NPRNA"])
        self.vmix = float(params["VMIX"])
        self.rwat = float(params["RWAT"])
        self.rmax = float(params["RMAX"])
        self.dx = float(params["DX"])
        self.dr1 = float(params["DR1"])
        self.test = float(params["TEST"])
        self.teste = float(params["TESTE"])
        self.testy = float(params["TESTY"])
        self.testv = float(params["TESTV"])

        self.atoms = [Atom.from_dict(at) for at in parse_atoms(atomstr, atomconfstr)]

    def dumps(self) -> str:
        if self.jobnam is None:
            raise ValueError("KGRN: 'jobnam' has to be given!")

        params = self.to_dict()
        atomstr = "\n".join(format_atom_line(at.to_dict()) for at in self.atoms)
        atomconf = "\n".join(format_atom_block(at.to_dict()) for at in self.atoms)

        if self._update_date:
            params["date"] = datetime.now()
        params["atoms"] = atomstr
        params["atomconf"] = atomconf
        s = TEMPLATE.format(**params)
        return s.replace(".0e", ".d")
