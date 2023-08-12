# emtolib

> Python tools for the EMTO package by L. Vitos et al.

This project is a collection of tools to make working with the EMTO package
easier. It was developed for **my personal** workflow, so it might not be
useful for you. Also check the official [pyemto project](https://github.com/hpleva/pyemto)
for a more general purpose tool.


## Installation

Run the following command to install `emtolib`:
```bash
python3 -m pip install git+ssh://git@github.com/dylanljones/emtolib.git
```


## CLI

There are a few commands that can be run from the command line:

```bash
Usage: emtolib [OPTIONS] COMMAND [ARGS]...

  emtolib 0.1.4 - Tools for the EMTO package by L. Vitos et al.

Options:
  --help  Show this message and exit.

Commands:
  atom        Atom configuration
  auxdirs     Create the auxillary directories in the given directories.
  checkdos    Checks the *.dos files in the given directories for
              unphysical...
  clear       Clears the output files in the given directories.
  conv        Greps for the convergence message in the *.prn files in the...
  diff        Get the difference between the *.dat files in the given...
  element     Get information about the given element.
  get         Gets the given value from the *.dat files in the given...
  grep        Greps for a pattern in the *.prn files in the given...
  hopfield    Extracts the Hopfield values from the *.prn files in the
              given...
  iter        Greps for the iteration number in the *.prn files in the
              given...
  running     Check how many slurm jobs are still running
  set         Sets the given value in the *.dat files in the given...
  set-header  Sets the header of the *.dat files in the given directories.
  set_paths   Sets the given value in the *.dat files in the given...
  submit      Batch-run the EMTO simulations in the given directories using...
```

#### Examples

As an example, we will use the following directory structure:

```
app/
├── Nb
│   ├── Nb1
│   │   ├── nb.dat
│   │   ├── run_emto
│   │   ├── ...
│   ├── Nb2
│   │   ├── ...
├── V
...
```

We can use the `conv` command to check for convergence in all subdirectories:

```bash
> emtolib conv app/Nb

app/Nb/Nb1:  Converged in 36 iterations at 13:23  09-Aug-23
app/Nb/Nb2:  Converged in 45 iterations at 13:28  09-Aug-23
...
```

## Usage

### Configuration

To use more advanced features of the library, you will need to create a
configuration file. This file should be named ``emto.ini`` and placed in the
root directory of your EMTO calculations.

Minimal example of the ``emto.ini`` configuration file:

```ini
[general]
# Put your custom configuration here...
app = app
...

[emto]
# Optional: custom EMTO root
root = /path/to/emto
# Paths are relative to the EMTO root
emto2 = kgrn2/kgrn_cpa
...

[slurm]
# Email is required for SLURM jobs
mail_user = <your email>
...
```
