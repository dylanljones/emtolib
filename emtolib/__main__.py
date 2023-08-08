# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-25
"""This module defines the command line interface for emtolib."""

import click
import subprocess
from pathlib import Path
from emtolib.directory import walk_emtodirs, diff_emtodirs
from emtolib.errors import DOSReadError
from emtolib.files import generate_makefile
from emtolib.common import elements, WorkingDir
from emtolib.config import CONFIG
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
            click.echo(frmt_file(folder.path))
            for line in lines:
                line = frmt_grep_line(line.strip(), pattern)
                click.echo(f"  {line}")


@click.group(
    name="emtolib",
    help=f"emtolib {__version__} - Tools for the EMTO package by L. Vitos et al.",
)
def cli():
    pass


@cli.command(help="Greps for a pattern in the *.prn files in the given directories.")
@click.argument("pattern")
@click.option("--last", "-l", is_flag=True, default=False)
@click.option("--first", "-f", is_flag=True, default=False)
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def grep(pattern, first, last, recursive, paths):
    _grep(pattern, first, last, recursive, paths)


@cli.command(
    name="iter",
    help="Greps for the iteration number in the *.prn files in the given directories.",
)
@click.option("--last", "-l", is_flag=True, default=False)
@click.option("--first", "-f", is_flag=True, default=False)
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def iter_command(first, last, recursive, paths):
    _grep("Iteration", first, last, recursive, paths)


@cli.command(help="Greps for Convergence in the *.prn files in the given directories.")
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def conv(recursive, paths):
    _grep("Converged", first=False, last=True, recursive=recursive, paths=paths)


@cli.command(help="Gets the given value from the *.dat files in the given directories.")
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("key", type=str, nargs=1)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def get(recursive, key, paths):
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        click.echo(f"{path} {key}={dat[key]}")


@cli.command(
    name="set", help="Sets the given value in the *.dat files in the given directories."
)
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("value", type=str, nargs=1)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def set_cmd(recursive, value, paths):
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    key, val = value.split("=")
    key, val = key.strip(), val.strip()
    for folder in folders:
        path = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        click.echo(f"{path} Setting {key} to {val}")
        dat[key] = val
        dat.dump()


@cli.command(
    help="Checks the *.dos files in the given directories for unphysical values."
)
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def checkdos(recursive, paths):
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


@cli.command(
    help="Generate a makefile for running all simulations in the given directory."
)
@click.argument("path", type=click.Path(), nargs=1, required=False)
def makefile(path):
    path = Path(path)
    click.echo(f"Generating makefile for directories in {path}")
    make = generate_makefile(path)
    make.dump()


@cli.command(
    help="Get the difference between the *.dat files in the given directories."
)
@click.option("--only_keys", "-k", is_flag=True, default=False)
@click.argument("path", type=click.Path(), nargs=1, required=False)
def diff(only_keys, path):
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


@cli.command(help="Clears the output files in the given directories.")
@click.option("--aux", "-a", is_flag=True, default=False)
@click.option("--keep", "-k", is_flag=True, default=False)
@click.argument("path", type=click.Path(), nargs=1, required=False)
def clear(aux, keep, path):
    folders = list(walk_emtodirs(path, recursive=True))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Clearing folder")
        folder.clear(aux=aux, keep=keep)


@cli.command(help="Sets the header of the *.dat files in the given directories.")
@click.option("--header", "-h", type=str, default="")
@click.option("--frmt", "-f", type=str, default="%d %b %y")
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def set_header(header, frmt, recursive, paths):
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        dat = folder.dat
        click.echo(f"{p} Setting header: '{header}'")
        dat.set_header(header, frmt)
        dat.dump()


@cli.command(help="Get information about the given element.")
@click.argument("symbol", type=click.Path(), nargs=1, required=True)
@click.argument("keys", type=str, nargs=-1, required=False)
def element(symbol, keys):
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


@cli.command(help="Create the auxillary directories in the given directories.")
@click.option("--keep", "-k", is_flag=True, default=False)
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def auxdirs(keep, recursive, paths):
    folders = list(walk_emtodirs(*paths, recursive=recursive))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        p = frmt_file(f"{str(folder.path) + ':':<{maxw}}")
        click.echo(f"{p} Creating aux directories")
        folder.mkdirs(keep=keep)


@cli.command(help="Batch-run the EMTO simulations in the given directories.")
@click.option("--executable", "-x", type=str, default="executable")
@click.option("--recursive", "-r", is_flag=True, default=False)
@click.argument("paths", type=click.Path(), nargs=-1, required=False)
def runbatched(executable, recursive, paths):
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
        slurm.set_body(executable, folder.dat.path.name)
        # Run slurm
        with WorkingDir(folder.path):
            cmd = f"sbatch {slurm.path.name}"
            stdout = subprocess.check_output(cmd, shell=True)
            stdout = stdout.decode("utf-8").replace("\n", "")
            click.echo(f"{p} {stdout}")


if __name__ == "__main__":
    cli()
