# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-21

import re
from ..common import EmtoFile
from ..config import update_slurm_settings

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
        mail_user=None,
        partition="epyc",
        ntasks=1,
        nodes=1,
        mail_type="FAIL,END,INVALID_DEPEND,TIME_LIMIT",
        time="7-0",
        mem="2gb",
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

    def update(self, **params):
        if "jobname" in params:
            self.jobname = params["jobname"]
        if "partition" in params:
            self.partition = params["partition"]
        if "ntasks" in params:
            self.ntasks = params["ntasks"]
        if "nodes" in params:
            self.nodes = params["nodes"]
        if "mail_type" in params:
            self.mail_type = params["mail_type"]
        if "mail_user" in params:
            self.mail_user = params["mail_user"]
        if "time" in params:
            self.time = params["time"]
        if "mem" in params:
            self.mem = params["mem"]

    def update_from_config(self, executable="", input_file="", conf=None):
        update_slurm_settings(self, executable, input_file, conf)
