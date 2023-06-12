# Python EMTO tools


## Installation

```bash
pip install git+ssh://git@github.com/dylanljones/emtolib.git
```


## Usage

### Configuration

Example ``emto.ini`` configuration file:

```ini
[emto]

root = ~/EMTO
kstr = EMTO5.8/kstr/smx
bmdl = EMTO5.8/bmdl/mdl
executable = EMTO5.8.1/kgrn/kgrn_cpa
executable2 = EMTO5.8.2/kgrn/kgrn_cpa
executable_dmft = EMTO_DMFT/kgrn/kgrn_cpa


[slurm]

partition = epyc
ntasks = 1
nodes = 1
mail_user = name@example.com
mail_type = FAIL,END,INVALID_DEPEND,TIME_LIMIT
mem = 2gb
time = 7-0
```
