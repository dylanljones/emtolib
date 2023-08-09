# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-25
"""This module defines the command line interface for emtolib."""

import click
import functools
import subprocess
from pathlib import Path
from emtolib.directory import walk_emtodirs, diff_emtodirs, EmtoDirectory
from emtolib.errors import DOSReadError
from emtolib.files import generate_makefile
from emtolib.common import elements, WorkingDir
from emtolib.config import CONFIG, update_emto_paths
from emtolib import __version__


def frmt_file(s):
    return click.style(s, fg="blue")


def frmt_header(s, maxw=0, color="blue"):
    return click.style(f"{str(s) + ':':<{maxw}}", fg=color)


def frmt_grep_line(line, pattern):
    line = line.strip()
    if pattern in line:
        line = line.replace(pattern, click.style(pattern, fg="green"))
    return line


def error(s):
    return click.style(s, fg="red")


def _grep(pattern, first, last, recursive, paths):
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        prn = folder.prn
        if not prn:
            continue
        lines = prn.grep(pattern).strip().splitlines()
        if not lines:
            continue
        if first or last:
            path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
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
    _grep(pattern, first, last, recursive, paths)


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
    _grep("Iteration", False, not all, recursive, paths)


@cli.command()
@multi_path_opts
def conv(recursive, paths):
    """Greps for the convergence message in the *.prn files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    pattern = "Converged"
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        prn = folder.prn
        if not prn:
            continue
        lines = prn.grep(pattern).strip().splitlines()
        if lines:
            line = frmt_grep_line(lines[-1].strip(), pattern)
            click.echo(f"{path} {line}")
        else:
            lines = prn.grep("Iteration").strip().splitlines()
            if not lines:
                click.echo(f"{path} {error('Not converged')}")
            else:
                line = lines[-1].strip()
                line = line.replace("KGRN:  Iteration no.", "")
                iteration = int(line.split()[0])
                click.echo(f"{path} {error('Not converged')} (Iter: {iteration})")


@cli.command()
@click.option("--dmft", "-d", is_flag=True, default=False, help="Use DMFT input files.")
@click.argument("key", type=str, nargs=1)
@multi_path_opts
def get(dmft, recursive, key, paths):
    """Gets the given value from the *.dat files in the given directories.

    KEY: The key of the value to get.
    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
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
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    key, val = value.split("=")
    key, val = key.strip(), val.strip()
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dmft if dmft else folder.dat
        click.echo(f"{path} Setting {key} to {val}")
        dat[key] = val
        dat.dump()


@cli.command(name="set_paths")
@click.option("--kstr", "-k", type=str, help="KSTR file name", default=None)
@click.option("--bmdl", "-b", type=str, help="BMDL file name", default=None)
@click.option("--kstr2", "-K", type=str, help="KSTR2 file name", default="")
@multi_path_opts
def set_paths(kstr, bmdl, kstr2, recursive, paths):
    """Sets the given value in the *.dat files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
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
@multi_path_opts
def checkdos(recursive, paths):
    """Checks the *.dos files in the given directories for unphysical values.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
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


@cli.command()
@single_path_opts
def makefile(path):
    """Generate a makefile for running all simulations in the given directory.

    PATH: The root path containing the EMTO directories.
    """
    path = Path(path)
    click.echo(f"Generating makefile for directories in {path}")
    make = generate_makefile(path)
    make.dump()


@cli.command()
@click.option(
    "--only_keys", "-k", is_flag=True, default=False, help="Only show key as output."
)
@single_path_opts
def diff(only_keys, path):
    """Get the difference between the *.dat files in the given directories.

    PATH: The root path containing the EMTO directories.
    """
    root = Path(path)
    diffs = diff_emtodirs(root)
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
        maxww = [max(len(str(val)) for val in diff.keys()) for diff in diffs.values()]
        maxw = max(maxww) + 1
        for p, d in diffs.items():
            click.echo(frmt_file(p))
            for key, val in d.items():
                click.echo(f"  {key + '=':<{maxw}} {val}")


@cli.command(name="clear")
@click.option("--aux", "-a", is_flag=True, default=False, help="Also clear aux files.")
@multi_path_opts
def clear_cmd(aux, recursive, paths):
    """Clears the output files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Clearing folder")
        folder.clear(aux=aux)


@cli.command()
@click.option("--header", "-h", type=str, default="", help="The header to set.")
@click.option(
    "--frmt", "-f", type=str, default="%d %b %y", help="The date format of the header."
)
@multi_path_opts
def set_header(header, frmt, recursive, paths):
    """Sets the header of the *.dat files in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        click.echo(f"{p} Setting header: '{header}'")
        dat.set_header(header, frmt)
        dat.dump()


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
@click.option("--keep", "-k", is_flag=True, default=False, help="Create keep files.")
@multi_path_opts
def auxdirs(keep, recursive, paths):
    """Create the auxillary directories in the given directories.

    PATHS: One or multiple paths to search for EMTO directories.
    """
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Creating aux directories")
        folder.mkdirs(keep=keep)


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

    folders = list(walk_emtodirs(*paths, recursive=recursive))
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
        # Run slurm
        with WorkingDir(folder.path):
            cmd = f"sbatch {slurm.path.name}"
            stdout = subprocess.check_output(cmd, shell=True)
            stdout = stdout.decode("utf-8").replace("\n", "")
            click.echo(f"{p} {stdout}")


@cli.group(name="atom", help="Atom configuration")
def atom_group():
    pass


@atom_group.command(name="add")
@click.option(
    "--clear", "-c", is_flag=True, default=False, help="Clear existing atoms first."
)
@click.argument("symbol", type=str, nargs=1, required=True)
@click.argument("kwargs", type=str, nargs=-1, required=False)
@single_path_opts
def add_atom(clear, symbol, kwargs, path):
    """Add a new atom to the input file in the given directory

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


if __name__ == "__main__":
    cli()
