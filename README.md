# Python EMTO tools


## Installation

```bash
python3 -m pip install git+ssh://git@github.com/dylanljones/emtolib.git
```

## CLI

There are a few commands that can be run from the command line:

```bash
usage: emtolib v0.1.0 [-h] {grep,conv,iter,set,get,check_dos,makefile,diff} ...

positional arguments:
  {grep,conv,iter,set,get,check_dos,makefile,diff}
    grep                Grep emto directories
    conv                Grep converged
    iter                Grep iter
    set                 Set variable of EMTO input file
    get                 Get variable of EMTO input file
    check_dos           Check DOS output
    makefile            Create makefile to run all EMTO folders
    diff                Diff multiple EMTO folders

options:
  -h, --help            show this help message and exit
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


- `grep`:
    Same as the grep command, but with automatic EMTO directory lookup. Example:

    ```bash
    emtolib grep "total energy" app/Nb
    ``` 

- `conv`:
    Check if a calculation has converged. Example:

    ```bash
    emtolib conv app/Nb
    ```

- `iter`:
    Check the iterations of one or more calculations.

    ```bash
    emtolib iter app/Nb
    ```

- `set`:
    Set paramters of the EMTO input file(s) for all EMTO directories in a 
    given root directory.
    ```bash
    emtolib set NKY=25 app/Nb
    ```
  
- `get`:
    Get paramters of the EMTO input file(s) for all EMTO directories in a 
    given root directory.
    ```bash
    emtolib get NKY app/Nb
    ```

- `check_dos`:
    Check the DOS of one or more calculations for unphysical values.

    ```bash
    emtolib check_dos app/Nb
    ```

- `makefile`:
    Create a makefile to run all EMTO folders in a given root directory.

    ```bash
    emtolib makefile app/Nb
    ```

- `diff`:
    Diff all EMTO folders in a given root directory.

    ```bash
    emtolib diff app/Nb
    ```

## Usage

### Configuration

To use more advanced features of the library, you will need to create a
configuration file. This file should be named ``emto.ini`` and placed in the
root directory of your EMTO calculations.

Minimal example of the ``emto.ini`` configuration file:

```ini
[emto]

root = ~/EMTO
# The following paths are relative to `root`
kstr = kstr/smx
bmdl = bmdl/mdl
executable = kgrn/kgrn_cpa
executable2 = kgrn/kgrn_cpa
executable_dmft = kgrn/kgrn_cpa


[slurm]

# Your slurm settings
partition = epyc
ntasks = 1
nodes = 1
mail_user = name@example.com
mail_type = FAIL,END,INVALID_DEPEND,TIME_LIMIT
mem = 2gb
time = 7-0
```
