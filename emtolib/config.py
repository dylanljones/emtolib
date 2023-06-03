# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from configparser import ConfigParser
from .input_files import EmtoKgrnFile


def read_config(file="emto.ini"):
    parser = ConfigParser()
    parser.read(file)
    config = dict()

    root = Path(parser["emto"]["root"])
    config["emto"] = dict(
        root=root,
        executable=str(root / parser["emto"]["executable"]).replace("\\", "/"),
        executable2=str(root / parser["emto"]["executable2"]).replace("\\", "/"),
        executable_dmft=str(root / parser["emto"]["executable_dmft"]).replace(
            "\\", "/"
        ),
        kstr=str(root / parser["emto"]["kstr"]).replace("\\", "/"),
        bmdl=str(root / parser["emto"]["bmdl"]).replace("\\", "/"),
    )
    config["slurm"] = dict(parser["slurm"])

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
