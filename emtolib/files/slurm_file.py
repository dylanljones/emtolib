# coding: utf-8
#
# This code is part of emtolib.
#
# Copyright (c) 2023, Dylan Jones

import re
from ..common import EmtoFile

RE_SBATCH = re.compile("#SBATCH --(.*?)=(.*?)$")

DEFAULT_COMMANDS = r"""
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=false
export OMP_STACKSIZE=256M
ulimit -Ss unlimited

echo - Host machine:\ \ \ $HOSTNAME
echo - I/O directory:\ \ $PWD
echo - SCF directory:\ \ $SCFDIR
echo - Calculation started:\ \ `date`
echo - PATH:\ \ $PATH
echo - E_RAND:\ \ $E_RAND

time {executable} < {file}

echo
echo - Calculation finished: `date`

"""


class SlurmScript(EmtoFile):
    def __init__(
        self,
        path="",
        jobname=None,
        partition=None,
        ntasks=None,
        nodes=None,
        mail_type=None,
        mail_user=None,
        time=None,
        mem=None,
        commands=None,
    ):
        super().__init__(path)
        self.jobname = jobname
        self.partition = partition
        self.ntasks = ntasks
        self.nodes = nodes
        self.mail_type = mail_type
        self.mail_user = mail_user
        self.time = time
        self.mem = mem
        self.commands = list()

        self.load(missing_ok=True)
        if commands:
            if isinstance(commands, str):
                commands = commands.splitlines(keepends=False)
            self.commands = list(commands)

    def iter_commands(self):
        return enumerate(list(self.commands))

    def find_command(self, pattern):
        regex = re.compile(pattern)
        for i, line in enumerate(self.commands):
            if regex.match(line):
                return i, line

    def set_body(self, executable, file):
        body = DEFAULT_COMMANDS.format(executable=executable, file=file)
        self.commands = body.splitlines(keepends=False)

    def dumps(self) -> str:
        if not self.jobname:
            raise ValueError("Required field 'job-name' not found!")
        if not self.mail_user:
            raise ValueError("Required field 'mail-user' not found!")

        lines = list()
        lines.append("#!/bin/bash")  # shebang
        lines.append(f"#SBATCH --job-name={self.jobname}")
        if self.partition:
            lines.append(f"#SBATCH --partition={self.partition}")
        if self.ntasks:
            lines.append(f"#SBATCH --ntasks={self.ntasks}")
        if self.nodes:
            lines.append(f"#SBATCH --nodes={self.nodes}")
        if self.mail_type:
            lines.append(f"#SBATCH --mail-type={self.mail_type}")
        if self.mail_user:
            lines.append(f"#SBATCH --mail-user={self.mail_user}")
        if self.time:
            lines.append(f"#SBATCH --time={self.time}")
        if self.mem:
            lines.append(f"#SBATCH --mem={self.mem}")
        # lines.append("")
        for cmd in self.commands:
            lines.append(cmd)
        return "\n".join(lines)

    def loads(self, data: str) -> None:
        lines = data.splitlines(keepends=False)
        if lines[0].startswith("#!"):
            lines.pop(0)  # shebang

        sbatch = dict()
        commands = list()
        for line in lines:
            line = line.strip()
            match = RE_SBATCH.match(line)
            if match:
                key = match.group(1)
                val = match.group(2)
                sbatch[key] = val
            elif not line.startswith("#"):
                commands.append(line)

        self.jobname = sbatch.get("job-name")
        self.partition = sbatch.get("partition")
        self.ntasks = sbatch.get("ntasks")
        self.nodes = sbatch.get("nodes")
        self.mail_type = sbatch.get("mail-type")
        self.mail_user = sbatch.get("mail-user")
        self.time = sbatch.get("time")
        self.mem = sbatch.get("mem")

        self.commands = commands
