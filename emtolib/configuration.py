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
    config = dict()
    emto = parser["emto"]
    root = Path(emto["root"])
    config["emto"] = dict(
        root=root,
        executable=str(root / emto["executable"]).replace("\\", "/"),
        executable1=str(root / emto["executable1"]).replace("\\", "/"),
        executable2=str(root / emto["executable2"]).replace("\\", "/"),
        executable_dmft=str(root / emto["executable_dmft"]).replace("\\", "/"),
        kstr=str(root / emto["kstr"]).replace("\\", "/"),
        bmdl=str(root / emto["bmdl"]).replace("\\", "/"),
    )
    config["slurm"] = dict(parser["slurm"])
    config.update(**parser["general"])

    return config


def update_slurm_settings(slurm, config, executable, jobname):
    slurm.ntasks = config["slurm"]["ntasks"]
    slurm.nodes = config["slurm"]["nodes"]
    slurm.mail_user = config["slurm"]["mail_user"]
    slurm.mail_type = config["slurm"]["mail_type"]
    slurm.mem = config["slurm"]["mem"]
    slurm.time = config["slurm"]["time"]

    executable = executable.replace("\\", "/")
    i, _ = slurm.find_command("time")
    slurm.commands[i] = f"time {executable} < {jobname}.dat"


def update_emto_paths(dat: EmtoKgrnFile, config, kstr, bmdl, kstr2=""):
    kstr_path = config["emto"]["kstr"] + "/" + kstr
    if not kstr2:
        kstr2 = kstr
    kstr2_path = config["emto"]["kstr"] + "/" + kstr2
    bmdl_path = config["emto"]["bmdl"] + "/" + bmdl

    dat.for001 = kstr_path
    dat.for001_2 = kstr2_path
    dat.for004 = bmdl_path


try:
    config = read_config()
except FileNotFoundError:
    config = dict()
