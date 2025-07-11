# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-25

"""Command line interface for emtolib."""

import functools
import subprocess
from pathlib import Path

import click
import numpy as np

from emtolib import __version__
from emtolib.common import WorkingDir, elements
from emtolib.config import CONFIG, update_emto_paths
from emtolib.errors import DOSReadError
from emtolib.directory import EmtoDirectory, diff_emtodirs, walk_emtodirs
from emtolib.easy_kgrn import run_emto_local, run_emto_slurm

ITER_TMPLT = "Iteration: {iter:>3} Etot: {etot:11.6f} erren: {erren:10.8f}"


def frmt_file(s):
    return click.style(str(s), fg="magenta")


def frmt_header(s, maxw=0, color="magenta"):
    return click.style(f"{str(s) + ':':<{maxw}}", fg=color)


def frmt_grep_line(line, pattern):
    line = line.strip()
    if pattern in line:
        line = line.replace(pattern, click.style(pattern, fg="green"))
    return line


def error(s):
    return click.style(s, fg="red")


def get_emtodirs(*paths, recursive=False):
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    # Recursively sort the folders by path paths
    if not folders:
        raise click.ClickException("No EMTO directories found.")
    maxparts = max([len(f.path.parts) for f in folders])
    for i in reversed(range(maxparts)):
        folders = sorted(folders, key=lambda f: f.path.parts[i])
    return folders


def single_path_opts(func):
    """Click argument decorator for commands accepting a single input path."""

    @click.argument("path", type=click.Path(), nargs=1, required=False, default=".")
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


def multi_path_opts(func):
    """Click argument decorator for commands accepting multiple input paths."""

    @click.option(
        "--recursive",
        "-r",
        is_flag=True,
        default=False,
        help="Recursively search for EMTO directories.",
    )
    @click.argument("paths", type=click.Path(), nargs=-1, required=False)
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@click.group(
    name="emtolib",
    help=f"emtolib {__version__} - Tools for the EMTO package by L. Vitos et al.",
)
def cli():
    pass


@cli.command(name="update")
def update():
    import os
    import sys

    cmd = fr"{sys.executable} -m pip install git+https://github.com/dylanljones/emtolib.git"

    click.echo(f"Updating emtolib: Running command '{cmd}'")
    click.echo()
    os.system(cmd)


# ======================================================================================
#                                   PRN files
# ======================================================================================


@cli.command(name="grep")
@click.argument("pattern")
@click.option("--last", "-l", is_flag=True, default=False, help="Only show last line")
@click.option("--first", "-f", is_flag=True, default=False, help="Only show first line")
@multi_path_opts
def grep_cmd(
    pattern,
    first,
    last,
    recursive,
    paths,
):
    """Greps for a pattern in the *.prn files in the given directories.

    PATTERN: The pattern to search for.
    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        prn = folder.prn
        if not prn:
            click.echo(f"{path} {error('No *.prn file found!')}")
            continue
        lines = prn.grep(pattern).strip().splitlines()
        if not lines:
            continue
        if first or last:
            if first:
                line = frmt_grep_line(lines[0].strip(), pattern)
                click.echo(f"{path} {line}")
            if last:
                line = frmt_grep_line(lines[-1].strip(), pattern)
                click.echo(f"{path} {line}")
        else:
            click.echo(frmt_file(str(folder.path)))
            for line in lines:
                line = frmt_grep_line(line.strip(), pattern)
                click.echo(f"  {line}")


@cli.command(name="iter")
@click.option("--all", "-a", is_flag=True, default=False, help="Show all line")
@multi_path_opts
def iter_command(
    all,  # noqa
    recursive,
    paths,
):
    """Greps for the iteration number in the *.prn files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        prn = folder.prn
        if not prn:
            click.echo(f"{path} {error('No *.prn file found!')}")
            continue
        iterations = prn.get_iterations()
        if not all:
            if not iterations:
                click.echo(f"{path} {error('No iterations!')}")
            else:
                it = iterations[-1]
                click.echo(f"{path} {ITER_TMPLT.format(**it)}")
        else:
            click.echo(frmt_file(str(folder.path)))
            if not iterations:
                click.echo(f"  {error('No iterations!')}")
            else:
                for it in iterations:
                    click.echo(f"  {ITER_TMPLT.format(**it)}")


@cli.command()
@multi_path_opts
def conv(recursive, paths):
    """Greps for the convergence message in the *.prn files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    pattern = "Converged"
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        prn = folder.prn
        if not prn:
            click.echo(f"{path} {error('No *.prn file found!')}")
            continue
        lines = prn.grep(pattern).strip().splitlines()
        if lines:
            line = frmt_grep_line(lines[-1].strip(), pattern)
            click.echo(f"{path} {line}")
        else:
            iterations = prn.get_iterations()
            if not iterations:
                click.echo(f"{path} {error('Not converged')}")
            else:
                it = iterations[-1]
                line = ITER_TMPLT.format(**it)
                click.echo(f"{path} {error('Not converged')} ({line})")


@cli.command()
@click.option("--mean", "-m", is_flag=True, default=False, help="Show sublattice Hopfields (mean)")
@click.option("--sum", "-s", is_flag=True, default=False, help="Sum over spins")
@multi_path_opts
def hopfield(mean, sum, recursive, paths):  # noqa
    """Extracts the Hopfield values from the *.prn files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        prn = folder.prn
        if prn is None:
            click.echo(f"{path} {error('No *.prn file found')}")
            continue
        if mean:
            if folder.dat is None:
                click.echo(f"{path} {error('No *.dat file found')}")
                continue
            try:
                eta = prn.get_sublat_hopfields(folder.dat)
            except AssertionError:
                click.echo(f"{path} {error('No Hopfields found')}")
                continue
            if sum:
                eta = eta.sum(axis=1)
            click.echo(f"{path} {eta}")
        else:
            click.echo(frmt_file(str(folder.path)))
            hopfields = prn.extract_hopfields()
            if not hopfields:
                click.echo(f"  {error('No Hopfields found')}")
                continue
            for atom, eta in hopfields:
                eta = np.array(eta)
                if sum:
                    eta = eta.sum(axis=0)
                click.echo(f"  {atom:<3}  {eta}")


# ======================================================================================
#                      DAT files / Directory structure
# ======================================================================================


@cli.command()
@click.option("--dmft", "-d", is_flag=True, default=False, help="Use DMFT input files.")
@click.argument("key", type=str, nargs=1)
@multi_path_opts
def get(dmft, recursive, key, paths):
    """Gets the given value from the *.dat files in the given directories.

    KEY: The key of the value to get.
    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dmft if dmft else folder.dat
        click.echo(f"{path} {key}={dat[key]}")


@cli.command(name="set")
@click.option("--dmft", "-d", is_flag=True, default=False, help="Use DMFT input files.")
@click.argument("value", type=str, nargs=1)
@multi_path_opts
def set_cmd(dmft, recursive, value, paths):
    """Sets the given value in the *.dat files in the given directories.

    VALUE: The key of the value to set. Must be in the form KEY=VALUE.
    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    key, val = value.split("=")
    key, val = key.strip(), val.strip()
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dmft if dmft else folder.dat
        click.echo(f"{path} Setting {key} to {val}")
        dat[key] = val
        dat.dump()


@cli.command(name="set-paths")
@click.option("--kstr", "-k", type=str, help="KSTR file name", default=None)
@click.option("--bmdl", "-b", type=str, help="BMDL file name", default=None)
@click.option("--kstr2", "-K", type=str, help="KSTR2 file name", default="")
@multi_path_opts
def set_paths(kstr, bmdl, kstr2, recursive, paths):
    """Sets the given value in the *.dat files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        if dat is None:
            click.echo(f"{path} {error('No dat file')}")
        click.echo(f"{path} Setting paths")
        _kstr = dat.for001 if kstr is None else kstr
        _bmdl = dat.for004 if bmdl is None else bmdl
        _kstr2 = dat.for001_2 if kstr2 is None else kstr2
        click.echo(f"KSTR={_kstr} BMDL={_bmdl} KSTR2={_kstr2}")
        update_emto_paths(dat, _kstr, _bmdl, _kstr2)
        dat.dump()


@cli.command()
@click.option("--only_keys", "-k", is_flag=True, default=False, help="Only show key as output.")
@multi_path_opts
def diff(only_keys, recursive, paths):
    """Get the difference between the *.dat files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    diffs = diff_emtodirs(*paths, recursive=recursive)
    if not diffs:
        click.echo("No differences found")
        return
    maxw = max(len(str(path)) for path in diffs.keys()) + 1
    if only_keys:
        for path, d in diffs.items():
            p = frmt_file(f"{str(path) + ':':<{maxw}}")
            vals = ", ".join(d.keys())
            click.echo(f"{p} {vals}")
    else:
        maxww = [max(len(str(val)) for val in d.keys()) for d in diffs.values()]
        maxw = max(maxww) + 1
        for p, d in diffs.items():
            click.echo(frmt_file(p))
            for key in sorted(list(d.keys())):
                click.echo(f"  {key + '=':<{maxw}} {d[key]}")


@cli.command("set-header")
@click.option("--header", "-h", type=str, default="", help="The header to set.")
@click.option("--frmt", "-f", type=str, default="%d %b %y", help="The date format of the header.")
@multi_path_opts
def set_header(header, frmt, recursive, paths):
    """Sets the header of the *.dat files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        click.echo(f"{p} Setting header: '{header}'")
        dat.set_header(header, frmt)
        dat.dump()


@cli.command(name="clear")
@click.option("--aux", "-a", is_flag=True, default=False, help="Also clear aux files.")
@multi_path_opts
def clear_cmd(aux, recursive, paths):
    """Clears the output files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Clearing folder")
        folder.clear(aux=aux)


@cli.command(name="clear-slurm")
@multi_path_opts
def clear_slurm(recursive, paths):
    """Clears the SLURM scripts in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Clearing SLURM scripts")
        folder.clear_slurm()


@cli.command()
@click.option("--keep", "-k", is_flag=True, default=False, help="Create keep files.")
@multi_path_opts
def auxdirs(keep, recursive, paths):
    """Create the auxillary directories in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Creating aux directories")
        folder.mkdirs(keep=keep)


# == Atom configuration ================================================================


@cli.group(name="atom", help="Atom configuration")
def atom_group():
    pass


@atom_group.command(name="add")
@click.option("--clear", "-c", is_flag=True, default=False, help="Clear existing atoms first.")
@click.argument("symbol", type=str, nargs=1, required=True)
@click.argument("kwargs", type=str, nargs=-1, required=False)
@single_path_opts
def add_atom(clear, symbol, kwargs, path):
    """Add a new atom to the input file in the given directory.

    SYMBOL: The symbol of the atom to add
    KWARGS: The keyword arguments to pass to the atom constructor. Must be given as
            'key=value' pairs.
    PATH: The path to the directory containing the input file
    """
    folder = EmtoDirectory(path, missing_ok=True)
    if not folder.exists():
        click.echo(error(f"Directory {path} not found"))
        return
    dat = folder.dat
    if not dat:
        click.echo(error("Input file not found"))
        return
    if clear:
        click.echo(frmt_header(str(folder.path)) + " Clearing existing atoms")
        dat.atoms.clear()

    if kwargs:
        kwargs_dict = dict()
        for kwarg in kwargs:
            key, value = kwarg.split("=")
            kwargs_dict[key] = value
    else:
        kwargs_dict = dict()
    vals = ",".join(f"{key}={value}" for key, value in kwargs_dict.items())
    click.echo(frmt_header(str(folder.path)) + f" Adding atom {symbol} ({vals})")
    dat.add_atom(symbol, **kwargs_dict)
    dat.update_mnta()
    dat.dump()


# ======================================================================================
#                               DOS files
# ======================================================================================


@cli.command()
@multi_path_opts
def checkdos(recursive, paths):
    """Checks the *.dos files in the given directories for unphysical values.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        try:
            dosfile = folder.dos
        except DOSReadError:
            click.echo(f"{path} " + error("Could not not read DOS file"))
            continue
        try:
            energy, dos = dosfile.get_total_dos()
            # Check for unphysical (negative) DOS values
            if (dos < 0).any():
                click.echo(f"{path} " + error("Negative DOS values found"))
                continue
            else:
                click.echo(f"{path} " + click.style("DOS ok", fg="green"))
        except AttributeError:
            click.echo(f"{path} " + error("Could not not read DOS file"))
            continue


def split_dosfile(folder: EmtoDirectory, with_sublattice: bool = False) -> None:
    dosfile = folder.dos
    kwargs = dict(index=False, sep=" ", header=False, float_format="%.4f")

    tdos_indices = list(set(dosfile.tdos.index))
    for spin in tdos_indices:
        df = dosfile.get_tdos(spin)
        file = f"tdos_{spin.lower()}.dos"
        df.to_csv(folder.path / file, **kwargs)

    pdos_indices_raw = set(dosfile.pdos.index)
    if not with_sublattice:
        pdos_indices_raw = set((atom, spin) for _, atom, spin in pdos_indices_raw)
        pdos_indices = list(pdos_indices_raw)
        for atom, spin in pdos_indices:
            df = dosfile.get_pdos(atom=atom, spin=spin)
            file = f"pdos_{atom}_{spin.lower()}.dos"
            df.to_csv(folder.path / file, **kwargs)
    else:
        pdos_indices = list(pdos_indices_raw)
        for sublatt, atom, spin in pdos_indices:
            df = dosfile.get_pdos(sublatt=sublatt, atom=atom, spin=spin)
            file = f"pdos_{sublatt}_{atom}_{spin.lower()}.dos"
            df.to_csv(folder.path / file, **kwargs)


@cli.command()
@multi_path_opts
@click.option("--sublattice", "-s", is_flag=True, default=False, help="Also split into sublattices.")
def split_dos(recursive, paths, sublattice: bool = False):
    """Split the DOS file into separate files for each spin, atom and optionally sublattice.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        try:
            dosfile = folder.dos
        except DOSReadError:
            click.echo(f"{path} " + error("Could not not read DOS file"))
            continue
        try:
            click.echo(f"{path} Splitting DOS file")
            split_dosfile(folder, with_sublattice=sublattice)
        except AttributeError:
            click.echo(f"{path} " + error("Could not not read DOS file"))
            continue


# ======================================================================================
#                                   Other
# ======================================================================================


@cli.group(name="slurm", help="SLURM setup")
def slurm_group():
    pass


@slurm_group.command(name="get")
@click.argument("key", type=str, nargs=1)
@multi_path_opts
def get_slurm_cmd(recursive, key, paths):
    """Gets the given value from the *.slurm files in the given directories.

    KEY: The key of the value to get.
    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    key = key.lower()
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        slurm = folder.slurm
        click.echo(f"{path} {key}={slurm[key]}")


@slurm_group.command(name="set")
@click.argument("value", type=str, nargs=1)
@multi_path_opts
def set_slurm_cmd(recursive, value, paths):
    """Sets the given value in the *.dat files in the given directories.

    VALUE: The key of the value to set. Must be in the form KEY=VALUE.
    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    key, val = value.split("=")
    key, val = key.strip().lower(), val.strip()
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        slurm = folder.slurm
        params = dat.to_dict()
        if key == "jobname":
            _val = val.format(**params)
        else:
            _val = val
        click.echo(f"{path} Setting {key} to {_val}")
        slurm[key] = _val
        slurm.dump()

# ======================================================================================
#                                   Other
# ======================================================================================


@cli.command()
@click.argument("symbol", type=str, nargs=1, required=True)
@click.argument("keys", type=str, nargs=-1, required=False)
def element(symbol, keys):
    """Get information about the given element.

    SYMBOL: The symbol of the element to get information about.
    KEYS: The keys to get information about. If not given, all keys will be shown.
    """
    try:
        el = elements[symbol]
    except KeyError:
        click.echo(error(f"Element {symbol} not found"))
        return

    click.echo(f"Element {el.symbol}:")
    if not keys:
        keys = list(el.keys())
    keys = sorted(list(keys))
    maxw = max(len(key) for key in keys) + 1
    for key in keys:
        click.echo("  " + frmt_header(key, maxw) + f" {el[key]}")


@cli.command()
@click.option("--executable", "-x", type=str, default="", help="The EMTO executable.")
@multi_path_opts
def submit(executable, recursive, paths):
    """Batch-run the EMTO simulations in the given directories using SLURM.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    emto_config = CONFIG["emto"]
    if executable in emto_config:
        root = Path(emto_config["root"])
        executable = root / emto_config[executable]

    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        slurm = folder.slurm
        if not slurm:
            click.echo(f"{p} {error('No slurm file found')}")
            continue
        # Update slurm body
        if executable:
            slurm.set_body(executable, folder.dat.path.name)
            slurm.dump()
        # Assure aux dirs exist before submitting job
        folder.mkdirs(keep=False)
        # Run slurm
        with WorkingDir(folder.path):
            cmd = f"sbatch {slurm.path.name}"
            stdout = subprocess.check_output(cmd, shell=True)
            stdout = stdout.decode("utf-8").replace("\n", "")
            click.echo(f"{p} {stdout}")


@cli.command()
@click.option("--local", "-l", is_flag=True, default=False, help="Run un-batched")
@multi_path_opts
def run(local, recursive, paths):
    """Run the EMTO simulations in the given directories using easy_kgrn.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        if local:
            click.echo(f"{p} Running locally")
            run_emto_local(folder.path)
        else:
            jobid = run_emto_slurm(folder.path)
            click.echo(f"{p} Submitted bacth job {jobid}")


@cli.command()
@click.option("--dry", "-d", is_flag=True, default=False, help="Run without canceling.")
@multi_path_opts
def cancel(dry, recursive, paths):
    folders = get_emtodirs(*paths, recursive=recursive)
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        # Find job IDs
        paths = folder.path.glob("slurm-*.out")
        job_ids = sorted([int(p.stem.replace("slurm-", "")) for p in paths])
        # Cancel latest job
        if len(job_ids) == 0:
            click.echo(f"{p} {error('No jobs found')}")
            continue

        last_id = job_ids[-1]
        with WorkingDir(folder.path):
            cmd = f"scancel {last_id}"
            if dry:
                click.echo(f"{p} {cmd}")
            else:
                try:
                    stdout = subprocess.check_output(cmd, shell=True)
                    stdout = stdout.decode("utf-8").replace("\n", "")
                    click.echo(f"{p} {stdout}")
                except subprocess.CalledProcessError:
                    click.echo(f"{p} {error('Job not found')}")


@cli.command()
@click.option("--all", "-a", is_flag=True, default=False, help="Show squeue output.")
def running(all):  # noqa
    """Check how many slurm jobs are still running."""
    cmd = r'squeue -u $(whoami) -o "%.8i %.10P %.14j %.8T %.10M %.8N"'
    stdout = subprocess.check_output(cmd, shell=True)
    stdout = stdout.decode("utf-8")
    lines = stdout.splitlines(keepends=False)
    nlines = len(lines) - 1
    click.echo(f"Running jobs: {nlines}")
    if all:
        click.echo("")
        for line in lines:
            click.echo(line)


if __name__ == "__main__":
    cli()
