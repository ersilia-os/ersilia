# Installation

**This is work in progress!** We are working day and night to have a first running version of the software before the end of 2020. Do not try to run this code yet.

## Python package

We recommend working inside a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) environment.
```bash
conda create -n ersilia python=3.7
conda activate ersilia
```
Then, simply install with pip.
```bash
# from pypi (feature not ready yet)
pip install ersilia
# or latest from github
pip install git+https://github.com/ersilia-os/ersilia.git
```
You are done!

```bash
ersilia --help
```

### Set up

The first time you run the `ersilia` command, additional dependencies will be installed and containers/images will be generated. So please be patient if it takes a while...

The following will happen automatically. For your information, this is what Ersilia is going to do:

1. Install conda and git. Most likely you have those installed in your computer.
2. Install RDKit and BioPython.
3. Create a `$HOME/eos` folder. This folder will be used to store models and meta-data.
4. Create a conda environment called `eosbase-<VERSION>`. This environment will be used as a base by other models. It is a bare minumum environment for executing the `ersilia` command.
5. Build a docker image named `ersiliaos/eos-server:<VERSION>`. This image will be used as a base for most of the docker images related to models.
