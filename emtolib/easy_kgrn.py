# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

import os
import subprocess
from configparser import ConfigParser
from emtolib import KgrnFile, EmtoDirectory


def parse_value(x):
    try:
        return int(x)
    except ValueError:
        pass
    try:
        return float(x)
    except ValueError:
        pass
    return x


def read_atoms(parser):
    atom_sections = [k for k in parser.sections() if k.startswith("atom:")]
    atom_sections.sort(key=lambda x: int(x.split(":")[1]))
    atoms = list()
    for sec in atom_sections:
        atom = {k: parse_value(v) for k, v in parser.items(sec)}
        if "u" in atom:
            atom["u"] = [float(x) for x in atom["u"].split(",")]
        if "j" in atom:
            atom["j"] = [float(x) for x in atom["j"].split(",")]
        atom["n"] = [int(x) for x in atom["n"].split()]
        atom["kappa"] = [int(x) for x in atom["kappa"].split()]
        atom["occup"] = [float(x) for x in atom["occup"].split()]
        atom["valen"] = [float(x) for x in atom["valen"].split()]
        atoms.append(atom)
    return atoms


def read_config_file(path):
    parser = ConfigParser()
    parser.read(path)

    params = {k: parse_value(v) for k, v in parser.items("general")}
    if "dmft" in parser.sections():
        dmft = {k: parse_value(v) for k, v in parser.items("dmft")}
        params.update(dmft)
        params["dmft"] = True
    else:
        params["dmft"] = False

    if "avg" in params:
        params["ga"] = True
    else:
        params["ga"] = False

    atom = {k: parse_value(v) for k, v in parser.items("atom")}
    params.update(atom)

    atoms = read_atoms(parser)

    run_config = dict(parser.items("run"))
    slurm = dict(parser.items("slurm"))

    return params, atoms, run_config, slurm


def convert_kgrn_input(params, atoms, kgrn_path):
    assert kgrn_path.endswith(".dat")

    file = KgrnFile(kgrn_path)
    if file.exists():
        file.atoms.clear()

    file.force_dmft(params.pop("dmft", False))
    file.force_ga(params.pop("ga", False))
    file.update(params)
    for atom in atoms:
        file.add_atom(atom.pop("symb"), **atom)
    if "mnta" not in params:
        file.update_mnta()
    return file


def shell_source(script):
    """Sometimes you want to emulate the action of "source" in bash,
    settings some environment variables. Here is a way to do it."""
    pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE, shell=True)
    output = pipe.communicate()[0].decode("utf-8")
    env = dict()
    for line in output.splitlines():
        try:
            key, value = line.split("=", 1)
            env[key] = value
        except ValueError:
            pass
    # env = dict((line.split("=", 1) for line in output.splitlines()))
    os.environ.update(env)


def run_emto_local(root_dir=".", input_file=".input.ini"):
    # cache cwd and cd into root_dir
    cwd = os.getcwd()
    os.chdir(root_dir)

    # Read input .ini file and prepare input .dat file and directory
    params, atoms, run_config, slurm_config = read_config_file(input_file)
    name = params["jobnam"] + ".dat"
    kgrn_file = convert_kgrn_input(params, atoms, name)
    kgrn_file.dump()
    folder = EmtoDirectory(".")
    folder.mkdirs()


    # Init shell
    if "init_script" in run_config:
        shell_source(run_config["init_script"])

    # Build and run command
    executable = run_config["executable"]
    cmd = f"{executable} < {folder.dat.path.name}"
    result = subprocess.run(cmd, shell=True, check=True, text=True)

    # Restore workding dir
    os.chdir(cwd)
    return result


def run_emto_slurm(root_dir=".", input_file=".input.ini"):
    # cache cwd and cd into root_dir
    cwd = os.getcwd()
    os.chdir(root_dir)

    # Read input .ini file and prepare input .dat file and directory
    params, atoms, run_config, slurm_config = read_config_file(input_file)
    name = params["jobnam"] + ".dat"
    kgrn_file = convert_kgrn_input(params, atoms, name)
    kgrn_file.dump()
    folder = EmtoDirectory(".")
    folder.mkdirs()

    executable = run_config["executable"]
    kgrn_in = folder.dat.path.name

    slurm = folder.get_slurm("run_emto")
    slurm.update(**slurm_config)
    slurm.set_body(executable, kgrn_in)
    slurm.dump()

    # Submit batch job
    stdout = subprocess.check_output(f"sbatch {slurm.path.name}", shell=True)
    text = stdout.decode("utf-8").replace("\n", "")
    jobid = int(text.split()[-1])
    print(f"Submitted batch job {jobid} ({slurm.jobname})")

    # Restore workding dir
    os.chdir(cwd)

    return jobid
