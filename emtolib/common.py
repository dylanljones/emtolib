# coding: utf-8
#
# This code is part of emto.
#
# Copyright (c) 2023, Dylan Jones

import os
import re
import json

RE_KEYVAL = re.compile(r"([a-zA-Z()0-9]+).*?=.*?([a-zA-Z0-9_\-.]+)")

with open(os.path.join(os.path.dirname(__file__), "elements.json"), "r") as _fp:
    elements = json.load(_fp)


def parse_params(data: str) -> dict:
    return dict(RE_KEYVAL.findall(data))


def parse_filepath(line: str) -> str:
    return line.split("=")[1].strip()


class EmtoFile:
    def __init__(self, path):
        self.path = path

    def loads(self, data: str) -> None:
        pass

    def dumps(self) -> str:
        pass

    def load(self, file: str = "") -> None:
        file = file or self.path
        with open(file, "r") as fp:
            self.loads(fp.read())

    def dump(self, file: str = "") -> None:
        file = file or self.path
        with open(file, "w") as fp:
            fp.write(self.dumps())

    def __repr__(self):
        return f"<{self.__class__.__name__}({self.path})>"
