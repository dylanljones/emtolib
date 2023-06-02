# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

import logging
import re
import json
from pathlib import Path

logger = logging.getLogger("emtolib")
sh = logging.StreamHandler()
sh.setLevel(logging.DEBUG)
logger.addHandler(sh)
fmt = logging.Formatter(
    "[%(asctime)s] %(levelname)-8s - %(name)-25s - %(message)s",
    datefmt="%H:%M:%S"
)
sh.setFormatter(fmt)
logger.setLevel(logging.INFO)

RE_KEYVAL = re.compile(r"([a-zA-Z()0-9]+).*?=.*?([a-zA-Z0-9_\-.]+)")

# Atomic masses in 1u ~ 1.660 10-27 kg
# Debye temperatures (in K)
# https://www.knowledgedoor.com/2/elements_handbook/debye_temperature.html
with open(Path(__file__).parent / "elements.json", "r") as _fp:
    elements = json.load(_fp)


def parse_params(data: str) -> dict:
    return dict(RE_KEYVAL.findall(data))


def parse_filepath(line: str) -> str:
    return line.split("=")[1].strip()


class EmtoFile:
    """Base class for EMTO input and output files."""

    def __init__(self, path):
        self.path = Path(path)

    def loads(self, data: str) -> None:
        pass

    def dumps(self) -> str:
        pass

    def load(self, file: str = ""):
        file = file or self.path
        with open(file, "r") as fp:
            self.loads(fp.read())
        return self

    def dump(self, file: str = "") -> None:
        file = file or self.path
        with open(file, "w") as fp:
            fp.write(self.dumps())

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.path})>"
