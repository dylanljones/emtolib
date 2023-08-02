# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-21

import logging
from pathlib import Path
from configparser import ConfigParser

logger = logging.getLogger(__name__)

# Locations to search for the config file (in order):
# Current working directory, Home directory, emtolib directory
LOCATIONS = [Path.cwd(), Path.home(), Path(__file__).parent.parent]


def _find_config_file(filename):
    logger.info("Looking for config file '%s'", filename)

    # Look for the config file in the supported locations
    for loc in LOCATIONS:
        logger.debug("Checking directory '%s'", loc)
        file = loc / filename
        if file.exists():
            logger.info("Found config file '%s'", file)
            return file
    logger.info("Config file '%s' not found!", filename)
    # Raise an error if none of the above worked
    raise FileNotFoundError(
        f"Could not find the config files {filename} in any of the directories: " +
        ", ".join(str(loc) for loc in LOCATIONS)
    )


def _unixpath(path):
    """Convert a path to a unix-style path, regardless of the current OS."""
    return str(path).replace("\\", "/")


def _load_config(file):
    parser = ConfigParser()
    parser.read(file)
    conf = dict()
    # Parse general section
    if "general" in parser:
        general = parser["general"]
        for key, val in general.items():
            try:
                conf[key] = Path(val)
            except TypeError:
                conf[key] = val

    # Parse emto section
    section = dict(parser["emto"])  # required section!
    root = Path(section.pop("root"))  # required key!
    # Ensure the following keys are present:
    emto = {
        "root": _unixpath(root),
        "kstr": _unixpath(root / section.pop("kstr", "kstr/smx")),
        "bmdl": _unixpath(root / section.pop("bmdl", "bmdl/mdl")),
        "executable": _unixpath(root / section.pop("executable", "kgrn/kgrn_cpa")),
    }
    # Other keys are optional, but valid
    for key, val in section.items():
        emto[key] = _unixpath(root / val)
    conf["emto"] = emto

    # Parse slurm section
    if "slurm" in parser:
        conf["slurm"] = dict(parser["slurm"])
    else:
        conf["slurm"] = dict()
    return conf


def read_config():
    # Find and read the default config file
    file = _find_config_file("emto.ini")
    conf = _load_config(file)
    return conf


def update_slurm_settings(slurm, conf, executable="", jobname=""):
    slurm.ntasks = conf["slurm"]["ntasks"]
    slurm.nodes = conf["slurm"]["nodes"]
    slurm.mail_user = conf["slurm"]["mail_user"]
    slurm.mail_type = conf["slurm"]["mail_type"]
    slurm.mem = conf["slurm"]["mem"]
    slurm.time = conf["slurm"]["time"]

    if executable:
        executable = executable.replace("\\", "/")
        i, _ = slurm.find_command("time")
        slurm.commands[i] = f"time {executable} < {jobname}.dat"


def update_emto_paths(dat, conf, kstr, bmdl, kstr2="", pot="pot/", chd="chd/", tmp=""):
    kstr_path = conf["emto"]["kstr"] + "/" + kstr

    bmdl_path = conf["emto"]["bmdl"] + "/" + bmdl
    dat.for001 = kstr_path
    if kstr2:
        kstr2_path = conf["emto"]["kstr"] + "/" + kstr2
        dat.for001_2 = kstr2_path
    dat.for004 = bmdl_path

    dat.dir002 = pot
    dat.dir003 = pot
    dat.dir006 = ""
    dat.dir009 = pot
    dat.dir010 = chd
    dat.dir011 = tmp


try:
    CONFIG = read_config()
except (FileNotFoundError, KeyError):
    CONFIG = dict()
