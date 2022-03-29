<div align="center">

# Welcome to the **Ersilia Model Hub**!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)



[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

</div>

![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)

### Table of Contents:
1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Features](https://github.com/ersilia-os/ersilia#features)
3. [Getting started](https://github.com/ersilia-os/ersilia#getting-started)
4. [Quick Start](https://github.com/ersilia-os/ersilia#quick-start)
5. [Usage](https://github.com/ersilia-os/ersilia#usage)
6. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
7. [Roadmap](https://github.com/ersilia-os/ersilia#roadmap)
8. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
9. [About us](https://github.com/ersilia-os/ersilia#about-us) 

# Project Description

The Ersilia Model Hub is the main project of the [Ersilia Open Source Initiative](https://ersilia.io). The Ersilia Model Hub's objective is to provide a platform for a user-friendly deployment of Artificial Intelligence and Machine Learning (AI/ML) models for drug discovery. The Ersilia Model Hub collects and deploys AI/ML models currently scattered in several scientific papers, code repositories, and other supplementary materials. The platform will enable scientists without coding expertise to easily browse through thousands of libraries relevant to their research, which will lead to more precise model predictions. 

Most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who do not have the necessary expertise. Pharmaceutical companies with the right technologies are reluctant to develop this area because they do not see a massive return on investment. This gap becomes even more prominent in low-income country institutions where access to bioinformatic facilities or data science experts is scarce.

Ersilia Model Hub aims to fulfill the ultimate goal of Ersilia, which is to lower the barrier to drug discovery, encouraging academic groups and companies to pursue the development of new medicines following the principles of Open Science. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.
* Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)

# Features

| Feature | Description |
| --- | --- |
| Data driven | Ersilia technology achieves state-of-the-art performance thanks to the integration of chemical, genomic and biomedical text data. Our AI tools are trained on millions of data points collected from the scientific literature and are available at no cost. |
| Bioactivity signatures | Our tools are designed to facilitate the use of AI/ML tools. Scientists can browse a collection of models, choose the ones relevant to their research interests and run predictions without writing a single line of code. |
| User friendly | Our tools are designed to facilitate the use of AI/ML tools. Scientists can browse a collection of models, choose the ones relevant to their research interests and run predictions without writing a single line of code. |
| Open source | All our assets are released under a permissive open source license. This means the scientific community can review, contribute and improve our code, resulting in tools validated more extensively than in the traditional peer-review system. |

# Getting Started

## Requirements

Please make sure you have the right installation of python and the other libraries and packages on your computer.
1. [Python 3.7 and above](https://www.python.org/)
2.  [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) / [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) / [Miniconda](https://docs.conda.io/en/latest/miniconda.html): the Ersilia tool will try to use it in order to handle dependencies safely and efficiently. 
3. [Docker](https://www.docker.com/): Please install Docker to ensure that all our AI/ML assets will work smoothly in your local device.
4. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [Git LFS](https://git-lfs.github.com/): All of our code is hosted in the Ersilia Github Organization profile


## Installation

### Install on Linux and MacOSX
**Open a terminal. The best is to set up a Conda environment.**

```
1 #create a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment
2 conda create -n ersilia [Python 3.7](https://www.python.org/)
3 #activate the environment
4 conda activate ersilia
```
**Then, simply install the Ersilia Python package.**
```
1 #clone from github
2 https://github.com/ersilia-os/ersilia.git
3 cd ersilia
4 #install with pip (use -e for developer mode)
5 pip install -e .
```
**You should be done! Quickly check that the CLI works on your terminal.**
```
1 #see ersilia CLI options
2 ersilia --help
```
## The Isaura data lake
We highly recommend installation of the[Python 3.7](https://github.com/ersilia-os/isaura) data lake. With Isaura, you will be able to cache your model predictions (i.e. store them in your local computer). Isaura is a relatively light Python package:
```
1 #clone from github
2 git clone https://github.com/ersilia-os/isaura.git
3 cd isaura
4 #install with pip (use -e for developer mode)
5 pip install -e
```
## Install on Windows
>We are not testing Windows installation consistently. If you encounter problems, please reach out to us at hello@ersilia.io.
We recommend that you install Ubuntu on your Windows machine. This can be now done very easily with WSL. You will need at least Windows 10.
Open a Power Shell with Admin permissions and type:
```
1 wsl --install
```
Then simply install the [Ubuntu](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) terminal on Windows.
Inside the Ubuntu terminal, you can now follow the installation instructions [above](https://github.com/ersilia-os/ersilia###Install-on-Linux-and-MacOSX).

## Command-line only installation snippets for Ubuntu
Below you can find a few snippets that can help you install dependencies in a Ubuntu terminal.

**<details><summary>CLICK TO VIEW INSTALLATION</summary>**
    
### The gcc compiler
**Probably you have the gcc compiler installed already. Just in case:**
```
1 sudo apt install build-essential
```
### Conda package manager
If you don't have the Conda package manager yet, we suggest you install Miniconda:
```
1 mkdir -p ~/miniconda3
2 wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
3 bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
4 rm -rf ~/miniconda3/miniconda.sh
5 ~/miniconda3/bin/conda init bash
6 ~/miniconda3/bin/conda init zsh
```
### GitHub CLI
Once Conda is installed (see above), you can use it to install the fantastic GitHub CLI:
```
1 conda install gh -c conda-forge
```
Use the GitHub CLI to login. This may be helpful if you have contributor permissions at Ersilia. Type:
```
1 gh auth login
```
And then follow the instructions.
### Git LFS
Git Large File Storage (LFS) can be installed from Conda as well:
```
1 conda install git-lfs -c conda-forge
```
Activate Git LFS:
```
1 git-lfs install
```
</details>

# Quick Start
## Command-line interface
We provide a command-line interface (CLI) to interact with the [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP). 

To check the available commands, simply type:
```
1 #list available commands
2 ersilia --help
```
## Browse the model catalog
You can explore our catalog of models. The following will return a list of models currently available in our remote repositories.
```
1 #catalog of models
2 ersilia catalog
```
## Fetch model and install it locally
The first step is to download the model to your local device and install it along with its dependencies. By default, a `~/eos` directory (for Ersilia Open Source) will be created in your `HOME`. This folder will contain all fetched models along with additional files to manage the AI/ML content available locally.

To download and install the model, simply use the `fetch` command:
```
1 #fetch model from remote repository
2 ersilia fetch <model>
```
## Serve model
Once the model has been fetched, it should be ready to be used. A model in the [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP) can be thought of as a set of APIs. You can serve the model like this:
```
1 ersilia serve <model>
```
# Usage
Let's consider the chemprop-antibiotic as an example. First you need to fetch your model with the Ersilia CLI as shown above:
## Make predictions
```
1 ersilia fetch chemprop-antibiotic
```
Generate a few (5) example molecules, to be used as input. Molecules are typically expressed in SMILES format.
```
1 ersilia example chemprop-antibiotic -n 5 -f my_molecules.csv
```
Then, serve your model:
```
1 ersilia serve chemprop-antibiotic
```
And run the prediction API:
```
1 ersilia api -i my_molecules.csv -o my_predictions.csv
```
Finally, close the service when you are done.
```
1 ersilia close
```
Please see the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailed explanations.

## Close model
Once you are done with predictions, it is advised to stop the model server:
```
1 #close model
2 ersilia close
```
## Delete model
If you are sure you don't want to use a model anymore, you may want to remove it from your computer. This includes deleting all model files and specific dependencies:
```
1 #delete model
2 ersilia delete chemprop-antibiotic
```

# Contribute

The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
- If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
- If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
- If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

# Roadmap
We are working to grow the Hub organically and responding to our users needs. Here a detail of the next features to come, stay tuned!

1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)

# License and Citation
This repository is open-sourced under the MIT License. Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.

# About Us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https://github.com/opencollective.com/ersilia).




