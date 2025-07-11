---
description: Documentation to run models on-premises
---

# Local inference

The Ersilia Model Hub is conveniently offered as a python package through PyPi and CondaForge, and each model is individually packaged as a Docker container.&#x20;

## Installation in Linux/MacOS

Ersilia is only maintained for Linux and Mac Operating Systems. If you work in Windows please use a [Windows Subsystem Linux](https://learn.microsoft.com/en-us/windows/wsl/install).

#### Prerequisites

* Python: we maintain Ersilia for Python 3.8 and above. Please make sure you have a compatible Python version installed on your computer. Visit the [official Python site](https://www.python.org) to learn more.
* Conda: ensure either [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Miniconda](https://docs.anaconda.com/miniconda/) are available in your system. This is the command to install it in **Ubuntu** (the command may be different if you do not use Ubuntu):

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

* Docker: Docker containers are an excellent way to share applications and simplify the management of system requirements and configurations. Please [install Docker](https://www.docker.com) to ensure that all our AI/ML assets will work smoothly in your local device.

#### Install from PyPi

```bash
# install ersilia from PyPi in an independent conda environment
conda create -n ersilia python=3.12
conda activate ersilia
pip install ersilia
```

#### Install from CondaForge

```bash
# install ersilia from CondaForge in an independent conda environment
conda create -n ersilia python=3.12
conda activate ersilia
conda install -c conda-forge ersilia
```

Once the Ersilia Model Hub is installed, test that it works by running the --help command:

```bash
ersilia --help
```

## Model Usage

You can explore the available models through [our website](https://ersilia.io/model-hub) or by running the following command in the CLI:

```bash
# display ready to use models with its eos identifier and title
ersilia catalog --hub --more
```

Each model is identified by:

* EOS-ID: `eos[1-9][a-z0-9]{3}`
* Slug: 1-3 word reference for the model
* Title: brief description of the model

Additionally, each model is associated to a `Task` e.g. Regression, an `Output`, e.g. Score and an `Output Shape` e.g Single, or List.&#x20;

Throughout this documentation, we will use the model `eos2r5a` (retrosynthetic-accessibility) as an example. This model has been incorporated from the paper _Retrosynthetic accessibility score (RAscore) – rapid machine learned synthesizability classification from AI driven retrosynthetic planning_ by [Thakkar et al, 2021](http://dx.doi.org/10.1039/D0SC05401A). The RA score is particularly useful to pre-screen large libraries of compounds, for example those produced by generative models. RA scores lie in the \[0,1] range; higher values indicate greater confidence that a compound is synthetically accessible.

To use a model, there are a few basic commands:

```
ersilia fetch eos2r5a
ersilia serve eos2r5a
ersilia run -i input.csv -o output.csv
ersilia close
```

The fetch command will download the model from DockerHub. Please make sure to have Docker active in your system before fetching a model. The serve command will bring the model alive anytime you want to use it. With the run command, you can make predictions by specifying the input and output files. Finally, close the model.

{% hint style="info" %}
If you serve a model that is not available locally, Ersilia will try to fetch it automatically, from DockerHub first and then defaulting to S3 if Docker is not active. Alternatively, models can be fetched from Github as well.&#x20;
{% endhint %}

### Input and output

The Ersilia Model Hub takes **chemical structures** as input, which should be specified as SMILES strings. To obtain the SMILES string of your compounds, you can use resources like [PubChem](https://pubchem.ncbi.nlm.nih.gov/).

Ersilia only accepts an input file in csv format, with one column and a header. Predictions are returned in tabular format as either .csv or .h5 depending on what is specified by the user:

{% code title="input.csv" %}
```bash
smiles
C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```
{% endcode %}

```bash
# predict using an input file and output file
ersilia run -i input.csv -o output.csv
```

### Other interesting commands

You can also get more information through the model card:

```bash
# display model card using ersilia identifier
ersilia catalog --card eos2r5a
```

### Delete model

If you are sure you don't want to use a model anymore, you may want to remove it from your computer. This includes deleting all model files and specific dependencies:

```bash
# delete model
ersilia delete retrosynthetic-accessibility
# or use the eos identifier
ersilia delete eos2r5a
```

## As a Python package

Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called `ErsiliaAPI`:

```python
# import main class
from ersilia.api.create_api import ErsiliaAPI
# instantiate the model(ex: Retrosynthetic Accessibility Score)
mdl = ErsiliaModel("eos2r5a")
```

Then, you can perform the same actions as in the CLI. To fetch and serve:

```python
# fetch model
mdl.fetch()

#serve model
mdl.serve()
```

To make **predictions** for Halicin and Ibuprofen:

```python
# Halicin and Ibuprofen
input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
]
# predict
mdl.run(input, output=None, batch_size=3)
# Pass the input for the model as a .csv file with only one column with header or a list of SMILES compounds
# Flexibility to specify file path for desired output file or automatically create an output file in the current directory
# Specify the batch size for generating model predictions. By default, Ersilia works with batch size of 100 inputs.

```

To **close** the model:

```python
# close model
mdl.close()
```

### Using the `with` statement

A more concise way to run prediction would be to use the `with` clause:

```python
# use with statement
# this allows for automatic closing of model
with mdl as model:
    model.info()
```
