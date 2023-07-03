# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from configparser import ConfigParser
from .files import EmtoKgrnFile


def read_config(file="emto.ini"):
    parser = ConfigParser()
    parser.read(file)
    conf = dict()
    emto = parser["emto"]
    root = Path(emto["root"])
    conf["emto"] = dict(
        root=root,
        executable=str(root / emto["executable"]).replace("\\", "/"),
        executable1=str(root / emto["executable1"]).replace("\\", "/"),
        executable2=str(root / emto["executable2"]).replace("\\", "/"),
        executable_dmft=str(root / emto["executable_dmft"]).replace("\\", "/"),
        kstr=str(root / emto["kstr"]).replace("\\", "/"),
        bmdl=str(root / emto["bmdl"]).replace("\\", "/"),
    )
    conf["slurm"] = dict(parser["slurm"])
    conf.update(**parser["general"])

    return conf


def update_slurm_settings(slurm, conf, executable, jobname):
    slurm.ntasks = conf["slurm"]["ntasks"]
    slurm.nodes = conf["slurm"]["nodes"]
    slurm.mail_user = conf["slurm"]["mail_user"]
    slurm.mail_type = conf["slurm"]["mail_type"]
    slurm.mem = conf["slurm"]["mem"]
    slurm.time = conf["slurm"]["time"]

    executable = executable.replace("\\", "/")
    i, _ = slurm.find_command("time")
    slurm.commands[i] = f"time {executable} < {jobname}.dat"


def update_emto_paths(dat: EmtoKgrnFile, conf, kstr, bmdl, kstr2=""):
    kstr_path = conf["emto"]["kstr"] + "/" + kstr
    if not kstr2:
        kstr2 = kstr
    kstr2_path = conf["emto"]["kstr"] + "/" + kstr2
    bmdl_path = conf["emto"]["bmdl"] + "/" + bmdl

    dat.for001 = kstr_path
    dat.for001_2 = kstr2_path
    dat.for004 = bmdl_path


try:
    CONFIG = read_config()
except (FileNotFoundError, KeyError):
    CONFIG = dict()
