# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

from pathlib import Path
from configparser import ConfigParser


def read_config(file="emto.ini"):
    parser = ConfigParser()
    parser.read(file)
    config = dict()

    root = Path(parser["emto"]["root"])
    config["emto"] = dict(
        root=root,
        executable=root / parser["emto"]["executable"],
        executable2=root / parser["emto"]["executable2"],
        executable_dmft=root / parser["emto"]["executable_dmft"],
        kstr=root / parser["emto"]["kstr"],
        bmdl=root / parser["emto"]["bmdl"],
    )
    config["slurm"] = dict(parser["slurm"])

    return config


def update_slurm_settings(slurm, config):
    slurm.ntasks = config["slurm"]["ntasks"]
    slurm.nodes = config["slurm"]["nodes"]
    slurm.mail_user = config["slurm"]["mail_user"]
    slurm.mail_type = config["slurm"]["mail_type"]
    slurm.mem = config["slurm"]["mem"]
    slurm.time = config["slurm"]["time"]
