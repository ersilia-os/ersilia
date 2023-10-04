---
description: Ersilia Model Hub Installation Instructions
---

# Installation

## 1. Install on Windows

{% hint style="warning" %}
We are not consistently testing Windows installation. If you encounter problems, please reach out to us at [hello@ersilia.io](mailto:hello@ersilia.io).
{% endhint %}

We recommend that you install Ubuntu on your Windows machine. This can now be done very easily with [WSL](https://docs.microsoft.com/en-us/windows/wsl/install). You will need at least Windows 10.

Open a PowerShell with Admin permissions and type:

```
wsl --install
```

Then simply [install the Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) on Windows. Inside the Ubuntu terminal, you can now follow the installation instructions for Linux below.

## 2. Install on Linux and MacOSX

### 2.1. Pre-requisites

There are a few third-party **pre-requisites**, and you need to have them installed on your computer before proceeding with the installation of the Ersilia tool. Don't worry; these dependencies are very popular, well-documented, and continuously maintained. Most probably, you already know them.

If you already meet some of the pre-requisites, you may skip the corresponding steps.

#### Pre-requisite 1: The GCC compiler

Probably, you already have the GCC compiler installed. This is the command to install it in **Ubuntu** (the command may be different if you do not use Ubuntu):

```
sudo apt install build-essential
```

#### Pre-requisite 2: Python and Conda

We have tested our tools on **Python 3.7** and above. Please make sure you have the correct Python installation on your computer. Visit the [official Python site](https://www.python.org) to learn more.

Although _not_ strictly necessary, we encourage the use of **Conda environments**. If Conda is available on the local device, the Ersilia tool will try to use it to handle dependencies safely and efficiently. Conda can be installed following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Many people choose to use [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installers.

If you don't have the Conda package manager yet, we suggest you install Miniconda. This is the command to install it in **Ubuntu** (the command may be different if you do not use Ubuntu):

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

#### Pre-requisite 3: Git and GitHub CLI

All of our code is hosted on the [Ersilia GitHub Organization](https://github.com/ersilia-os) profile. The Ersilia tool needs to fetch code from GitHub, so make sure Git is installed on your device. You can follow [these instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

In addition, we highly recommend the fantastic [GitHub CLI](https://cli.github.com/manual/installation), which Ersilia uses if available.

Once Conda is installed (see above), you can use it to install GitHub CLI

```
conda install gh -c conda-forge
```

Use the GitHub CLI to log in. This may be helpful if you have contributor permissions at Ersilia. Type:

```
gh auth login
```

And then follow the instructions.

#### Pre-requisite 4: Git LFS

Some pre-trained models, especially the ones that contain many parameters, require [**Git Large File Storage** (LFS)](https://git-lfs.github.com/).

Git LFS can be installed from Conda:

```
conda install git-lfs -c conda-forge
```

Activate Git LFS:

```
git-lfs install
```

#### Pre-requisite 5: The Isaura Data Lake

We highly recommend the installation of the [Isaura](https://github.com/ersilia-os/isaura) data lake. With Isaura, you will be able to cache your model predictions (i.e., store them on your local computer). Isaura is a relatively lightweight Python package.

Before installing Isaura, please ensure that you have created the Ersilia Conda environment. If you haven't done so, please follow the [instructions](#22-install-ersilia) below.

Once you have the Ersilia Conda environment ready, proceed to install Isaura:

```bash
# activate Ersilia's Conda environment
conda activate ersilia
python -m pip install isaura==0.1
```

{% hint style="danger" %}
Isaura has different functionalities; please make sure to install version 0.1 (v0.1).
{% endhint %}

#### Pre-requisite 6: Docker

Docker containers are an excellent way to share applications and simplify the management of system requirements and configurations. Please [install Docker](https://www.docker.com) to ensure that all our AI/ML assets will work smoothly on your local device.

### 2.2. Install Ersilia

After all the pre-requisites are met, we are ready to install the Ersilia tool.

Open a terminal and set up a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) environment.

```bash
# create a Conda environment
conda create -n ersilia python=3.10
# activate the environment
conda activate ersilia
```

Then, simply install the Ersilia Python package.

```bash
# clone from GitHub
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
# install with pip (use -e for developer mode)
pip install -e .
```

You should be done! Check that the CLI works on your terminal, and explore the available commands

```bash
# see Ersilia CLI options
ersilia --help

# see Ersilia's model catalog
ersilia catalog
```

{% hint style="info" %}
The Ersilia Model Hub is continually growing to fulfil the needs of the community. Please do not hesitate to [request new models](https://www.ersilia.io/request-model)!
{% endhint %}

###