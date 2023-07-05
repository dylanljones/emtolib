# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

import shutil
import re
import json
import logging
from pathlib import Path
from typing import Union
from collections import abc
import numpy as np

logger = logging.getLogger("emtolib")
sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)
logger.addHandler(sh)
fmt = logging.Formatter(
    "[%(asctime)s] %(levelname)-8s - %(name)-25s - %(message)s", datefmt="%H:%M:%S"
)
sh.setFormatter(fmt)
logger.setLevel(logging.INFO)

RE_KEYVAL = re.compile(r"([a-zA-Z()0-9]+).*?=.*?([a-zA-Z0-9_\-.]+)")


class Mapping(abc.Mapping):
    def __init__(self, data):
        self._data = data

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __getitem__(self, item):
        return self._data[item]

    def __getattr__(self, item):
        return self.__getitem__(item)

    def __repr__(self):
        return repr(self._data)

    def __str__(self):
        return str(self._data)


class ElementView(Mapping):
    def __repr__(self):
        return f"<ElementView {self.symbol} ({self.name})>"

    def __str__(self):
        data = self._data.copy()
        name = data.pop("name")
        symbol = data.pop("symbol")
        lines = list()
        lines.append(f"{symbol} ({name})")
        for key, value in data.items():
            lines.append(f"   {key + ':':<15} {value}")
        return "\n".join(lines) + "\n"


class Elements(Mapping):
    def __getitem__(self, item):
        data = self._data[item]
        return ElementView(data)

    def __repr__(self):
        return f"<Elements ({len(self)})>"


# Atomic masses in 1u ~ 1.660 10-27 kg
# Debye temperatures (in K)
# https://www.knowledgedoor.com/2/elements_handbook/debye_temperature.html
with open(Path(__file__).parent / "elements.json", "r") as _fp:
    elements = Elements(json.load(_fp))


def parse_params(data: str) -> dict:
    return dict(RE_KEYVAL.findall(data))


def parse_filepath(line: str) -> str:
    return line.split("=")[1].strip()


class PostInitCaller(type):
    """Metaclass to call __post_init__ after __init__."""

    def __call__(cls, *args, **kwargs):
        obj = type.__call__(cls, *args, **kwargs)
        obj.__post_init__()
        return obj


class EmtoFile(metaclass=PostInitCaller):
    """Base class for EMTO input and output files."""

    def __init__(self, path: Union[str, Path] = None):
        if path is None:
            path = ""
        self.path = Path(path)

    def __post_init__(self):
        """Called after __init__."""
        pass

    def loads(self, data: str) -> None:
        pass

    def dumps(self) -> str:
        pass

    def load(self, file: Union[str, Path] = "", missing_ok: bool = False):
        """Load data from file."""
        file = Path(file or self.path)
        if file.is_dir():
            if missing_ok:
                return None
            raise IsADirectoryError(f"{file} is a directory!")
        try:
            with open(file, "r") as fp:
                data = fp.read()
        except FileNotFoundError:
            if missing_ok:
                return None
            raise
        self.loads(data)
        return self

    def dump(self, file: Union[str, Path] = "") -> None:
        """Dump data to file."""
        file = file or self.path
        # Encode contents to UTF-8 and replace DOS with UNIX line endings
        # This is probably not necessary, but it's done anyway for good measure
        data = self.dumps().encode("utf-8").replace(b"\r\n", b"\n")
        with open(file, "wb") as fp:
            fp.write(data)

    def exists(self) -> bool:
        """Check if file exists."""
        return self.path.exists()

    def move(self, dst: Union[str, Path], exist_ok: bool = False) -> None:
        """Move file to new location."""
        dst = Path(dst)
        if dst.exists() and not exist_ok:
            raise FileExistsError(f"Destination {dst} already exists!")
        shutil.move(self.path, dst)
        self.path = dst

    def rename(self, name: str, exist_ok: bool = False) -> None:
        """Rename file."""
        self.move(self.path.parent / name, exist_ok)

    def copy(self, dst: Union[str, Path], exist_ok: bool = False):
        """Copy file to new location."""
        dst = Path(dst)
        if dst.exists():
            if not exist_ok:
                raise FileExistsError(f"Destination {dst} already exists!")
        else:
            shutil.copytree(self.path, dst)
        return self.__class__(dst)

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.path})>"


def dict_diff(d1, d2, exclude=None):
    # copy to avoid changing the original dict
    d1, d2 = d1.copy(), d2.copy()
    # Ensure exclude is a tuple
    if exclude is not None and isinstance(exclude, str):
        exclude = (exclude,)
    # Convert lists and arrays inside the dicts to tuples
    for k, v in d1.items():
        if isinstance(v, list) or isinstance(v, np.ndarray):
            d1[k] = tuple(v)
    for k, v in d2.items():
        if isinstance(v, list) or isinstance(v, np.ndarray):
            d2[k] = tuple(v)
    # Find the difference between the two dicts
    diffset = set(d1.items()) ^ set(d2.items())
    diff = dict()
    for key in dict(diffset).keys():
        if exclude is None or key not in exclude:
            diff[key] = (d1[key], d2[key])
    return diff
