# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-03

from pathlib import Path  # noqa: F401
from .common import logger, elements
from .config import read_config, CONFIG, update_emto_paths, update_slurm_settings
from .files import (
    KgrnFile,
    Atom,
    BmdlFile,
    PrnFile,
    DosFile,
    read_dos,
    load_dos,
    SlurmScript,
)
from .directory import EmtoDirectory, is_emtodir, walk_emtodirs
from .errors import (
    EmtoError,
    KGRNError,
    KGRNReadError,
    KGRNWriteError,
    DMFTError,
    DMFTReadError,
    DMFTWriteError,
    DOSError,
    DOSReadError,
)

try:
    from ._version import version as __version__
except ImportError:
    __version__ = "unknown"
