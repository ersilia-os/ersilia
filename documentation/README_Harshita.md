<div align="center">

   <a href="https://ersilia.io/" target="_blank"><img src="https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png" alt="Ersilia logo"></a>
    <br />

# A charity developing open source tools to facilitate global health drug discovery, with a focus on neglected diseases.


<a href="https://github.com/ersilia-os/ersilia/network/members"><img src="https://img.shields.io/github/forks/ersilia-os/ersilia" alt="Forks Badge"/></a>
<a href="https://github.com/ersilia-os/ersilia/graphs/contributors"><img alt="GitHub contributors" src="https://img.shields.io/github/contributors/ersilia-os/ersilia?color=2b9348"></a>
<a href="https://github.com/ersilia-os/ersilia"><img src="https://badges.frapsoft.com/os/v2/open-source.svg" alt="Open Source"/></a>
<a href="https://github.com/ersilia-os/ersilia/issues"><img src="https://img.shields.io/github/issues/ersilia-os/ersilia" alt="Issues Open"/></a>
[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012)[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Build Status](https://img.shields.io/github/workflow/status/ersilia-os/ersilia/Tests?label=tests&style=flat-square)](https://github.com/ersilia-os/ersilia/actions)
[![Twitter Account](https://img.shields.io/twitter/follow/ersiliaio?color=00acee&label=twitter&style=flat-square)](https://twitter.com/ersiliaio)
[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

</div>
<br>

## Now the big question is , What is Ersilla Model Hub 🤔 and how is it realated to The Ersilia Open Source ? 

Let's firstly understand about **Ersilla Model Hub**, So Ersilia Model Hub is a free, online, open-source platform where scientists can browse through a catalogue of
AI/ML models, select the ones that are relevant to their research and run online predictions without the need to write a single line of code. We are gathering, 
in a single resource, two classes of models.

***It is one of the largest repository of AI/ML models for Global Health.Our platform is aimed at helping researchers identify
drug candidates for orphan and neglected diseases, design molecules de novo, understand mechanisms of action or anticipate adverse side effects.***
<br> 

![image](https://user-images.githubusercontent.com/79797000/162241886-6d3d73bf-b892-4820-b4f3-f2f9250295a3.png)

<br>

Talking about the **Ersilia Open Source Initiative** , It is a UK-based charity whose mission is to strengthen research capacity against infectious and neglected diseases.
Currently under which there are total of 6 projects.The projects are as follows :

| Project          | Description |
| ----------------- | ------------------------------------------------------------------ |
| Capacity building | Resources and courses tailor-made to the needs of our collaborators |
| Automated modelling for low-resourced settings | End-to-end package to train and deploy AI/ML models based on chemistry data.|
| Ersilia Model Hub | An open-source repository of AI/ML models for global health and neglected diseases.|
| Virtual screening of African natural products | Analysis of the natural product landscape to identify new anti infectives.|
| Privacy-preserving AI | Development of an encryption tool to enable AI/ML model training using IP-sensitive data.|
| Open Source Malaria Series 4 | Design of new antimalarials using reinforcement learning algorithms.|

Find out more at: [https://ersilia.io/](https://ersilia.io/)

<details>
<summary>Table of Contents</summary>

- [Getting Started](#getting-started)
  - [Prerequisites](#Prerequisites)
  - [Installation](#installation)
    - [Linux and MacOSX](#linux-and-macosx)
    - [The Isaura data lake](#the-isaura-data-lake)
    - [Windows](#windows)
    - [Ubuntu](#ubuntu)
    - [Gcc compiler](#gcc-Compiler)
    - [Conda package manager](#conda-package-manager)
    - [GitHub CLI](#gitHub-CLI)
    - [Git LFS](#gitHub-LFS)


</details>

## Overview 📑

We collect models developed by third parties and available in scientific publications. On the other hand, we develop models in-house and/or in 
collaboration with research groups in LMICs. In other words, part of our philanthropic mission is to increase visibility and facilitate access to AI/ML research 
developed by the community, and part is to contribute AI/ML tools ourselves in order to fulfill unmet global health needs.  The Ersilia Model Hub currently works 
on a CLI interface

![image](https://user-images.githubusercontent.com/79797000/162250521-919c17ef-78a7-4f6d-bd30-dc17bbac2942.png)

<p align="right">(<a href="#top">back to top</a>)</p>

## Getting Started

### Prerequisites

Let us peep and set up the Prerequisites 🔖

1.	[x] [Python 3.7 and above](https://www.python.org/)
2.	[x] [Docker](https://www.docker.com/)
3.	[x] [Conda](https://www.docker.com/)   ***Anaconda or Miniconda can also be installed***
4.	[x] [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [Git LFS](https://git-lfs.github.com/)
5.	[x] A Unix-like operating system: macOS, Linux

After completing the Prerequisites let's move towards our next step ; Installation

### Installation 

<hr>


#### Linux and MacOSX

<hr>

**Open a terminal. The best is to set up a Conda environment.**

- Create a conda environment
```
conda create -n ersilia python=3.7
```

- Activate the environment
```
conda activate ersilia
```
**Then, simply install the Ersilia Python package.**

- Clone from github

```
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
```

- Install with pip (use -e for developer mode)

```
pip install -e .
```
**You should be done! Quickly check that the CLI works on your terminal.**

- See ersilia CLI options
```
ersilia --help
```
<hr>

#### The Isaura data lake
We highly recommend installation of the [Isaura](https://github.com/ersilia-os/isaura) data lake. With Isaura, you will be able to cache your model predictions (i.e. store them in your local computer). Isaura is a relatively light Python package:

- Clone from github
```
git clone https://github.com/ersilia-os/isaura.git
cd isaura
```
- Install with pip (use **-e** for developer mode)
```
pip install -e .
```
<hr>

#### Windows

*We recommend that you install Ubuntu on your Windows machine. This can be now done very easily with WSL. You will need at least Windows 10.*

- Open a Power Shell with **Admin permissions** and type:
```
wsl --install
```
Then simply [install the Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) on Windows.
Inside the Ubuntu terminal, you can now follow the installation instructions above.
<hr>

## Ubuntu

<hr>

Below you can find a few snippets that can help you install dependencies in a Ubuntu terminal.

<hr>
#### Gcc compiler

Probably you have the gcc compiler installed already. Just in case:

```
sudo apt install build-essential
```
<hr>

#### Conda package manager
If you don't have the Conda package manager yet, we suggest you install Miniconda:

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```
<hr>

#### GitHub CLI
Once Conda is installed (see above), you can use it to install the fantastic GitHub CLI:

```
conda install gh -c conda-forge
```

Use the GitHub CLI to login. This may be helpful if you have contributor permissions at Ersilia. Type:

```
gh auth login
```

And then follow the instructions.
<hr>

#### Git LFS
- Git Large File Storage (LFS) can be installed from Conda as well:

```
conda install git-lfs -c conda-forge
```
- Activate Git LFS:
```
git-lfs install
```

**You can also follow the **installation instructions** from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation).**

#### Model Generation

Once Ersilia is installed, you can **browse models** in the [Ersilia Model Hub](https://airtable.com/shrXfZ8pqro0jjcsG/tblZGe2a2XeBxrEHP/viwd5XJVLslkE11Tg).

Select one model. For example `chemprop-antibiotic`. You can **fetch** your model with the Ersilia CLI:
```
ersilia fetch chemprop-antibiotic
```
Generate a few (5) example molecules, to be used as input. Molecules are typically expressed in SMILES format.
```
ersilia example chemprop-antibiotic -n 5 -f my_molecules.csv
```
Then, **serve** your model:
```
ersilia serve chemprop-antibiotic
```
And run the prediction **API**:
```
ersilia api -i my_molecules.csv -o my_predictions.csv
```
Finally, **close** the service when you are done.
```
ersilia close
```

Please see the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailed explanations.
<p align="right"><a href="#top">Back to top</a></p>








## Our Vision 
A world with egalitarian access to healthcare


