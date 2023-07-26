# -*- coding: utf-8 -*
# Author: Dylan Jones
# Date:   2023-07-25

import sys
from argparse import ArgumentParser
from .files import KGRNError
from .directory import walk_emtodirs, is_emtodir, EmtoDirectory


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


def add_path_arg(parser):
    helpstr = "Path to one or more EMTO directories"
    parser.add_argument("paths", type=str, nargs="*", default=[""], help=helpstr)


def init_argparser():
    parser = ArgumentParser("emtolib")
    subparsers = parser.add_subparsers(dest="command")

    # Grep command
    parser_grep = subparsers.add_parser("grep", help="Grep emto directories")
    parser_grep.add_argument("pattern", type=str, help="Pattern to grep for")
    parser_grep.add_argument("-t", "--type", type=str, help="File types")
    add_path_arg(parser_grep)

    # Converged command
    parser_conv = subparsers.add_parser("converged", help="Grep converged")
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

    return parser


def handle_grep(args):
    for folder in iter_emtodirs(args):
        prn = folder.prn
        line = prn.grep(args.pattern).strip()
        if line:
            print(f"{folder.path}:\n{line.strip()}")


def handle_converged(args):
    for folder in iter_emtodirs(args):
        prn = folder.prn
        line = prn.grep("Converged", ignore_case=False).strip()
        if line:
            print(f"{folder.path}: {line}")
        else:
            print(f"{folder.path}: Not converged")


def handle_iter(args):
    for folder in iter_emtodirs(args):
        prn = folder.prn
        line = prn.grep("Iter", ignore_case=False).strip()
        if line:
            print(f"{folder.path}: {line}")


def handle_set(args):
    argvalue = args.value[0]
    key, value = argvalue.split("=")
    key = key.strip()
    value = value.strip()
    for folder in iter_emtodirs(args):
        dat = folder.dat
        print(f"{folder.path}: Setting {key} to {value}")
        dat[key] = value
        dat.dump()


def handle_get(args):
    key = args.key[0]
    for folder in iter_emtodirs(args):
        dat = folder.dat
        print(f"{folder.path}: {dat[key]}")


def handle_check_dos(args):
    for folder in iter_emtodirs(args):
        dosfile = folder.dos
        try:
            energy, dos = dosfile.get_total_dos()
            # Check for unphysical (negative) DOS values
            if (dos < 0).any():
                print(f"{folder.path}: Negative DOS values found")
                continue
            else:
                print(f"{folder.path}: DOS ok")
        except AttributeError:
            print(f"{folder.path}: No DOS file found")
            continue


HANDLERS = {
    "grep": handle_grep,
    "converged": handle_converged,
    "iter": handle_iter,
    "set": handle_set,
    "get": handle_get,
    "check_dos": handle_check_dos,
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
