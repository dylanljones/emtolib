# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-17

from ..common import EmtoFile, fort2py_float, py2fort_float
from ..ftmplt import Template


TEMPLATE = """\
#LMTO structure constants, number Matsubaras, Temp (Kelvin), Solver, DC, Mixing Sigma, Nomis
{for001:}
{nmats:d}
{temp:g}
{solver}
{dc:d}
{mixsig:g}
{nomis:d}
"""  # noqa


class DmftFile(EmtoFile):

    extension = ".dat"
    template = Template(TEMPLATE, ignore_case=True)

    def __init__(self, path, **kwargs):
        super().__init__(path)
        self.for001 = "bcc.dat"
        self.nmats = 1024
        self.temp = 400.0
        self.solver = "uppsalasolver"
        self.dc = 1
        self.mixsig = 0.1
        self.nomis = 80

        self.load(missing_ok=True)
        if kwargs:
            self.update(kwargs)

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
        return data

    def update(self, *args, **kwargs):
        data = dict(*args, **kwargs)
        for k, v in data.items():
            if not hasattr(self, k):
                raise KeyError(f"{k} is not a valid field of {self.__class__.__name__}")
            if isinstance(v, str):
                v = v.strip()
            self.__setattr__(k, v)

    def loads(self, text):
        text = fort2py_float(text)
        data = self.template.parse(text)
        self.update(data)
        return self

    def dumps(self) -> str:
        params = self.to_dict()
        return py2fort_float(self.template.format(params))
