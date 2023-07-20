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
    general = parser["general"]
    conf["data_dir"] = Path(general["data_dir"])
    conf["xarr_dir"] = Path(general["xarr_dir"])
    conf["fig_dir"] = Path(general["fig_dir"])
    conf["exp_dir"] = Path(general["exp_dir"])
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
