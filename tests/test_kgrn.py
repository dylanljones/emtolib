# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-03

from pathlib import Path
from pytest import mark
from numpy.testing import assert_array_equal
from emtolib.files import KgrnFile

TEST_ROOT = Path(__file__).parent.parent / ".testdata"

NB_AT_PARAMS = {
    "symbol": "Nb",
    "norb": 14,
    "ion": 0,
    "config": "4d3_5s2",
    "iq": 1,
    "it": 1,
    "ita": 1,
    "nz": 41,
    "conc": 1.0,
    "sms": 1.0,
    "sws": 1.0,
    "wswst": 1.0,
    "qtr": 0.0,
    "splt": 0.0,
    "fix": "N",
    "u": [],
    "j": [],
}

V_AT_PARAMS = {
    "symbol": "V",
    "norb": 9,
    "ion": 0,
    "config": "3d3_4s2",
    "iq": 1,
    "it": 1,
    "ita": 2,
    "nz": 23,
    "conc": 1.0,
    "sms": 1.0,
    "sws": 1.0,
    "wswst": 1.0,
    "qtr": 0.0,
    "splt": 0.0,
    "fix": "N",
    "u": [],
    "j": [],
}


def test_empty_file():
    _ = KgrnFile()
    _ = KgrnFile("test.dat")


@mark.parametrize(
    "path",
    [
        TEST_ROOT / "CPA" / "Nb" / "nb.dat",
        TEST_ROOT / "CPA" / "Nb25" / "nb.dat",
        TEST_ROOT / "CPA" / "Fe040_DLM" / "fe0484.dat",
        TEST_ROOT / "CPA" / "Sn50" / "co3sn2s2.dat",
    ],
)
def test_parse_format_cpa(path):
    text = path.read_text()
    dat = KgrnFile(path)
    text2 = dat.dumps()
    text = "\n".join(text.splitlines()[1:])
    text2 = "\n".join(text2.splitlines()[1:])
    assert text.strip() == text2.strip()


def test_parse_format_dmft_pure():
    path = TEST_ROOT / "DMFT" / "Nb" / "nb.dat"
    text = path.read_text()
    dat = KgrnFile(path)
    text2 = dat.dumps()
    text = "\n".join(text.splitlines()[1:])
    text2 = "\n".join(text2.splitlines()[1:])
    assert text.strip() == text2.strip()


def test_parse_format_dmft():
    path = TEST_ROOT / "DMFT" / "Nb25" / "nb.dat"
    text = path.read_text()
    dat = KgrnFile(path)
    text2 = dat.dumps()
    text = "\n".join(text.splitlines()[1:])
    text2 = "\n".join(text2.splitlines()[1:])
    assert text.strip() == text2.strip()


def test_parse_1():
    params = {
        "jobnam": "nb",
        "header": "03 Jun 23",
        "comment": "Self-consistent KKR calculation for nb",
        "strt": "N",
        "msgl": 1,
        "expan": "S",
        "fcd": "Y",
        "func": "SCA",
        "niter": 999,
        "nlin": 31,
        "nprn": 0,
        "ncpa": 7,
        "nt": 1,
        "mnta": 1,
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
        "nky": 57,
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
        "imagz": 0.001,
        "eps": 0.2,
        "elim": -1.0,
        "amix": 0.01,
        "efmix": 1.0,
        "vmtz": 0.0,
        "mmom": 0.0,
        "tole": 1e-07,
        "tolef": 1e-07,
        "tolcpa": 1e-06,
        "tfermi": 500.0,
        "sws": 3.071,
        "nsws": 1,
        "dsws": 0.05,
        "alpcpa": 0.602,
        "efgs": 0.1,
        "hx": 0.3,
        "nx": 7,
        "nz0": 16,
        "stmp": "N",
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
        "for001": "~/EMTO/EMTO5.8/kstr/smx/bcc.tfh",
        "for001_2": "~/EMTO/EMTO5.8/kstr/smx/bcc.tfh",
        "for004": "~/EMTO/EMTO5.8/bmdl/mdl/bcc.mdl",
        "dir002": "pot/",
        "dir003": "pot/",
        "dir006": "",
        "dir009": "pot/",
        "dir010": "chd/",
        "dir011": "/tmp/",
    }

    path = TEST_ROOT / "CPA" / "Nb" / "nb.dat"
    dat = KgrnFile(path)
    for k, v in params.items():
        assert dat[k] == v
        assert getattr(dat, k) == v

    assert params == dat.to_dict()

    assert len(dat.atoms) == 1
    at = dat.atoms[0]
    for k, v in NB_AT_PARAMS.items():
        assert getattr(at, k) == v

    assert_array_equal(at.n, [1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5])
    assert_array_equal(at.kappa, [-1, -1, 1, -2, -1, 1, -2, 2, -3, -1, 1, -2, 2, -1])
    assert_array_equal(at.occup, [2, 2, 2, 4, 2, 2, 4, 4, 6, 2, 2, 4, 3, 2])
    assert_array_equal(at.valen, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1])


def test_parse_2():
    expected = {
        "jobnam": "nb",
        "header": "03 Jun 23",
        "comment": "Self-consistent KKR calculation for nb",
        "strt": "N",
        "msgl": 1,
        "expan": "S",
        "fcd": "Y",
        "func": "SCA",
        "niter": 200,
        "nlin": 31,
        "nprn": 0,
        "ncpa": 7,
        "nt": 1,
        "mnta": 2,
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
        "nky": 57,
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
        "imagz": 0.001,
        "eps": 0.2,
        "elim": -1.0,
        "amix": 0.01,
        "efmix": 1.0,
        "vmtz": 0.0,
        "mmom": 0.0,
        "tole": 1e-07,
        "tolef": 1e-07,
        "tolcpa": 1e-06,
        "tfermi": 500.0,
        "sws": 2.8775,
        "nsws": 1,
        "dsws": 0.05,
        "alpcpa": 0.602,
        "efgs": 0.1,
        "hx": 0.3,
        "nx": 7,
        "nz0": 16,
        "stmp": "N",
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
        "for001": "~/EMTO/EMTO5.8/kstr/smx/bcc.tfh",
        "for001_2": "~/EMTO/EMTO5.8/kstr/smx/bcc.tfh",
        "for004": "~/EMTO/EMTO5.8/bmdl/mdl/bcc.mdl",
        "dir002": "pot/",
        "dir003": "pot/",
        "dir006": "",
        "dir009": "pot/",
        "dir010": "chd/",
        "dir011": "/tmp/",
    }

    path = TEST_ROOT / "CPA" / "Nb25" / "nb.dat"
    dat = KgrnFile(path)
    for k, v in expected.items():
        assert dat[k] == v
        assert getattr(dat, k) == v

    assert expected == dat.to_dict()

    assert len(dat.atoms) == 2
    at = dat.get_atom("Nb")
    assert at.symbol == "Nb"
    for k, v in NB_AT_PARAMS.items():
        if k == "conc":
            v = 0.25
        assert getattr(at, k) == v
    assert_array_equal(at.n, [1, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5])
    assert_array_equal(at.kappa, [-1, -1, 1, -2, -1, 1, -2, 2, -3, -1, 1, -2, 2, -1])
    assert_array_equal(at.occup, [2, 2, 2, 4, 2, 2, 4, 4, 6, 2, 2, 4, 3, 2])
    assert_array_equal(at.valen, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1])

    at = dat.get_atom("V")
    assert at.symbol == "V"
    for k, v in V_AT_PARAMS.items():
        if k == "conc":
            v = 0.75
        assert getattr(at, k) == v
    assert_array_equal(at.n, [1, 2, 2, 2, 3, 3, 3, 3, 4])
    assert_array_equal(at.kappa, [-1, -1, 1, -2, -1, 1, -2, 2, -1])
    assert_array_equal(at.occup, [2, 2, 2, 4, 2, 2, 4, 3, 2])
    assert_array_equal(at.valen, [0, 0, 0, 0, 0, 0, 0, 1, 1])
