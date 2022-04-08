<div id="top" align="center">

# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](https://github.com/ersilia-os/ersilia/blob/master/CODE_OF_CONDUCT.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

</div>

![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)

### Table of Contents:
1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Getting started](https://github.com/ersilia-os/ersilia#getting-started)
3. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
4. [Roadmap](https://github.com/ersilia-os/ersilia#roadmap)
5. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
6. [About us](https://github.com/ersilia-os/ersilia#about-us)
# Project Description
The Ersilia Model Hub is the main project of the [Ersilia Open Source Initiative](https://ersilia.io). The aim is to provide a platform for a user-friendly deployment of AI/ML models, where scientists can browse through the assets, identify those which is relevant to their research and obtain predictions without the need of coding expertise. Currently, most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who does not have the necessary expertise. This gap becomes even larger in Low and Middle Income Country institutions where access to bioinformatic facilities or data science experts are scarce. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.

# Overview

![image-1](https://user-images.githubusercontent.com/63330165/160850205-9d269457-06ad-46b7-9aaa-7a934c2fb47c.png)

- **Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)**</br>
- **Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)**

# Getting started

## Requirements

Please make sure you have the installation of below libraries on your machine.
1. [Python 3.7 and above](https://www.python.org/)
2. [Docker](https://www.docker.com/)
3.  [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) / [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) / [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
4. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [Git LFS](https://git-lfs.github.com/)

## Installation

### Installation on Linux and MacOSX

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
### The Isaura data lake
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

### Installation on Windows
*We recommend that you install Ubuntu on your Windows machine. This can be now done very easily with WSL. You will need at least Windows 10.*

- Open a Power Shell with **Admin permissions** and type:
```
wsl --install
```
Then simply [install the Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) on Windows.
Inside the Ubuntu terminal, you can now follow the installation instructions above.

## Command-line only installation snippets for Ubuntu

Below you can find a few snippets that can help you install dependencies in a Ubuntu terminal.
## The gcc compiler

Probably you have the gcc compiler installed already. Just in case:

```
sudo apt install build-essential
```
## Conda package manager
If you don't have the Conda package manager yet, we suggest you install Miniconda:

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

## GitHub CLI
Once Conda is installed (see above), you can use it to install the fantastic GitHub CLI:

```
conda install gh -c conda-forge
```

Use the GitHub CLI to login. This may be helpful if you have contributor permissions at Ersilia. Type:

```
gh auth login
```

And then follow the instructions.
## Git LFS
- Git Large File Storage (LFS) can be installed from Conda as well:

```
conda install git-lfs -c conda-forge
```
- Activate Git LFS:
```
git-lfs install
```

**You can also follow the **installation instructions** from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation).**

# Model Generation

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

# Features

| Feature | Description |
| --- | --- |
| **Bioactivity signatures** | Tools are designed to facilitate the use of AI/ML tools. Anyone even from non-coding background can browse a collection of models, choose the ones relevant to their research interests and run predictions without writing even a single line of code. |
| **Data driven** | Ersilia technology achieves state-of-the-art performance thanks to the integration of chemical, genomic and biomedical text data. |
| **User friendly platform** | Our tools are designed to facilitate the use of AI/ML tools and platform is user friendly, it is beginners-friendly as well. |
| **Open source** | All our assets are released under the MIT license. This means the anyone in ersilia community can review, contribute and improve our code, resulting in tools validated more extensively than in the traditional peer-review system. |

# Contribute
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
- If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
- If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
- If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

# Roadmap
We are working to grow the Hub organically and responding to our users needs. Here are the details of the next features to come, stay tuned!
1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)

# License and citation
This repository is open-sourced under the GNU AFFERO GENERAL PUBLIC LICENSE.
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.

# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https:/opencollective.com/ersilia).
<p align="right"><a href="#top">Back to top</a></p>
