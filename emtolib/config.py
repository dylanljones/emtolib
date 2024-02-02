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

# Default config values
CONFIG = {
    # EMTO paths: Should all be unix paths relative to the EMTO root
    "emto": {
        "root": "~/EMTO",
        "kstr": "kstr/smx",
        "bmdl": "bmdl/mdl",
        "emto": "kgrn/kgrn_cpa",
        "dmft": "kgrn_dmft/kgrn_cpa",
    },
    # General Slurm settings (without user-mail)
    "slurm": {
        "partition": "epyc",
        "ntasks": 1,
        "nodes": 1,
        "mail_type": "FAIL,END,INVALID_DEPEND,TIME_LIMIT",
        "mem": "2gb",
        "time": "7-0",
    },
}


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
        f"Could not find the config files {filename} in any of the directories: "
        ", ".join(str(loc) for loc in LOCATIONS)
    )


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
    conf["emto"] = dict(parser["emto"])
    # Parse slurm section
    conf["slurm"] = dict(parser["slurm"])

    return conf


def read_config(filename="emto.ini"):
    """Find and read the default config file."""
    file = _find_config_file(filename)
    conf = _load_config(file)
    return conf


def update_config(filename="emto.ini"):
    """Update the default config with the values from the config file."""
    logger.debug("Updating emtolib config")

    file = _find_config_file(filename)
    conf = _load_config(file)
    for key, values in conf.items():
        if isinstance(values, dict):
            if key in CONFIG:
                CONFIG[key].update(values)
            else:
                CONFIG[key] = dict(values)
        else:
            CONFIG[key] = values  # noqa


try:
    update_config()
except (FileNotFoundError, KeyError):
    pass


def update_emto_paths(
    dat, kstr, bmdl, kstr2="", pot=None, chd=None, tmp=None, conf=None
):
    if pot is None:
        pot = "pot/"
    if chd is None:
        chd = "chd/"
    if tmp is None:
        tmp = ""
    if not conf:
        conf = CONFIG
    emto_conf = conf["emto"]
    root = emto_conf["root"]
    kstr_path = root + "/" + emto_conf["kstr"] + "/" + kstr
    if kstr2:
        kstr2_path = root + "/" + emto_conf["kstr"] + "/" + kstr2
    else:
        kstr2_path = ""
    bmdl_path = root + "/" + emto_conf["bmdl"] + "/" + bmdl
    dat.update_paths(kstr_path, bmdl_path, kstr2_path, pot, chd, tmp)


def update_slurm_settings(slurm, executable="", input_file="", conf=None):
    if not conf:
        conf = CONFIG
    emto_conf = conf["emto"]
    root = emto_conf["root"]

    slurm.ntasks = conf["slurm"]["ntasks"]
    slurm.nodes = conf["slurm"]["nodes"]
    slurm.mail_user = conf["slurm"]["mail_user"]
    slurm.mail_type = conf["slurm"]["mail_type"]
    slurm.mem = conf["slurm"]["mem"]
    slurm.time = conf["slurm"]["time"]

    if executable:
        ex = emto_conf[executable]
        # executable = executable.replace("\\", "/")
        slurm.set_body(root + "/" + ex, input_file)
