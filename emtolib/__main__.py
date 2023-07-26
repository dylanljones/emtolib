# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-25

import sys
from argparse import ArgumentParser
from pathlib import Path
from .files import generate_makefile
from .errors import KGRNError
from .directory import walk_emtodirs, is_emtodir, EmtoDirectory, diff_emtodirs
from . import __version__


def iter_emtodirs(args):
    try:
        recursive = args.recursive
    except AttributeError:
        recursive = True
    for argpath in args.paths:
        try:
            if is_emtodir(argpath):
                yield EmtoDirectory(argpath)
        except KGRNError as e:
            print(e)

        for folder in walk_emtodirs(argpath, recursive=recursive):
            yield folder


def add_path_arg(parser, nargs="*", default=""):
    helpstr = "Path to one or more EMTO directories"
    parser.add_argument("paths", type=str, nargs=nargs, default=[default], help=helpstr)


def init_argparser():
    parser = ArgumentParser(f"emtolib v{__version__}")
    subparsers = parser.add_subparsers(dest="command")

    # Grep command
    parser_grep = subparsers.add_parser("grep", help="Grep emto directories")
    parser_grep.add_argument("pattern", type=str, help="Pattern to grep for")
    parser_grep.add_argument("-t", "--type", type=str, help="File types")
    add_path_arg(parser_grep)

    # Converged command
    parser_conv = subparsers.add_parser("conv", help="Grep converged")
    add_path_arg(parser_conv)

    # Iter command
    parser_iter = subparsers.add_parser("iter", help="Grep iter")
    add_path_arg(parser_iter)

    # Set command
    parser_set = subparsers.add_parser("set", help="Set variable of EMTO input file")
    parser_set.add_argument("value", type=str, nargs=1, help="Value to set (key=value)")
    add_path_arg(parser_set)

    # Set command
    parser_get = subparsers.add_parser("get", help="Get variable of EMTO input file")
    parser_get.add_argument("key", type=str, nargs=1, help="Key to get")
    add_path_arg(parser_get)

    # Check DOS command
    parser_set = subparsers.add_parser("check_dos", help="Check DOS output")
    add_path_arg(parser_set)

    # Makefile command
    parser_make = subparsers.add_parser(
        "makefile", help="Create makefile to run all EMTO folders"
    )
    add_path_arg(parser_make)

    # Diff command
    parser_diff = subparsers.add_parser("diff", help="Diff multiple EMTO folders")
    parser_diff.add_argument(
        "-x", "--exclude", nargs="*", type=str, help="Exclude keys"
    )
    parser_diff.add_argument("-k", "--only_keys", action="store_true", default=False)
    add_path_arg(parser_diff)

    # parser_set.add_argument(
    #     "-l", "--local", action="store_true", help="Local run (dont use slurm!)"
    # )

    return parser


def handle_grep(args):
    folders = list(iter_emtodirs(args))
    for folder in folders:
        path = folder.path
        prn = folder.prn
        line = prn.grep(args.pattern).strip()
        if line:
            print(f"{path}\n{line.strip()}")


def handle_converged(args):
    folders = list(iter_emtodirs(args))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = f"{str(folder.path) + ':':<{maxw}}"
        prn = folder.prn
        line = prn.grep("Converged", ignore_case=False).strip()
        if line:
            print(f"{path} {line}")
        else:
            print(f"{path} Not converged")


def handle_iter(args):
    folders = list(iter_emtodirs(args))
    for folder in folders:
        prn = folder.prn
        line = prn.grep("Iter", ignore_case=False).strip()
        if line:
            print(f"{folder.path}\n {line}")


def handle_set(args):
    folders = list(iter_emtodirs(args))
    maxw = max(len(str(folder.path)) for folder in folders) + 1

    argvalue = args.value[0]
    key, value = argvalue.split("=")
    key = key.strip()
    value = value.strip()
    for folder in folders:
        path = f"{str(folder.path) + ':':<{maxw}}"
        dat = folder.dat
        print(f"{path} Setting {key} to {value}")
        dat[key] = value
        dat.dump()


def handle_get(args):
    key = args.key[0]
    folders = list(iter_emtodirs(args))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = f"{str(folder.path) + ':':<{maxw}}"
        dat = folder.dat
        print(f"{path} {dat[key]}")


def handle_check_dos(args):
    folders = list(iter_emtodirs(args))
    maxw = max(len(str(folder.path)) for folder in folders) + 1
    for folder in folders:
        path = f"{str(folder.path) + ':':<{maxw}}"
        dosfile = folder.dos
        try:
            energy, dos = dosfile.get_total_dos()
            # Check for unphysical (negative) DOS values
            if (dos < 0).any():
                print(f"{path} Negative DOS values found")
                continue
            else:
                print(f"{path} DOS ok")
        except AttributeError:
            print(f"{path} No DOS file found")
            continue


def handle_makefile(args):
    if len(args.paths) > 1:
        raise ValueError("Only one root path allowed")

    path = Path(args.paths[0])
    make = generate_makefile(path)
    make.dump()


def handle_diff(args):
    if len(args.paths) > 1:
        raise ValueError("Only one root path allowed")

    root = Path(args.paths[0])
    diffs = diff_emtodirs(root, exclude=args.exclude)
    if not diffs:
        print("No differences found")
        return

    maxw = max(len(str(path)) for path in diffs.keys()) + 1
    if args.only_keys:
        for path, diff in diffs.items():
            p = f"{str(path) + ':':<{maxw}}"
            vals = ", ".join(diff.keys())
            print(f"{p} {vals}")
    else:
        maxww = [max(len(str(val)) for val in diff.keys()) for diff in diffs.values()]
        maxw = max(maxww) + 1
        for path, diff in diffs.items():
            p = str(path)
            vals = "\n".join(f"  {key + '=':<{maxw}} {val}" for key, val in diff.items())
            print(f"{p}\n{vals}")


HANDLERS = {
    "grep": handle_grep,
    "conv": handle_converged,
    "iter": handle_iter,
    "set": handle_set,
    "get": handle_get,
    "check_dos": handle_check_dos,
    "makefile": handle_makefile,
    "diff": handle_diff
}


def cli(argv):
    # Initialize argument parser and parse arguments
    parser = init_argparser()
    if not argv:
        parser.print_help(sys.stderr)
        return

    args = parser.parse_args(argv)
    try:
        handler = HANDLERS[args.command]
    except KeyError:
        raise ValueError(f"No handler found, unknown command: {args.command}")

    handler(args)


if __name__ == "__main__":
    cli(sys.argv[1:])
