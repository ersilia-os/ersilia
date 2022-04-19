<div id="top"></div>

# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)

### Table of Contents:
1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Getting started](https://github.com/ersilia-os/ersilia#getting-started)
3. [Operation](https://github.com/ersilia-os/ersilia#operation)
4. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
5. [Roadmap](https://github.com/ersilia-os/ersilia#roadmap)
6. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
7. [About us](https://github.com/ersilia-os/ersilia#about-us)

# Project Description
The Ersilia Model Hub is the main project of the [Ersilia Open Source Initiative] (https://ersilia.io). The aim is to provide a platform for a user-friendly deployment of AI/ML models, where scientists can browse through the assets, identify those which are relevant to their research and obtain predictions without the need of coding expertise. Currently, most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who do not have the necessary expertise. This gap becomes even larger in Low and Middle Income Country institutions where access to bioinformatic facilities or data science experts are scarce. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgment) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.

* Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)
<p align="right">(<a href="#top">back to top</a>)</p>

# Getting started
** Requirements**
*Python 3.7 and above
*MacOS or Linux or Windows 10

**Installation Instructions**
*For Linux and MacOS
- Set up Conda environment
```
conda create -n ersilia python=3.7

```
- Activate the environment

```
conda activate ersilia

```
- Clone the repository

```
git clone https://github.com/ersilia-os/ersilia.git

```
- Change the working directory
```
cd ersilia

```
- Install with pip
```
pip install -e

```
- Confirm CLI works on your terminal
```
ersilia --help

```

* Install Dependency [Isaura data lake is highly recommended]

- Clone from Github
```
git clone https://github.com/ersilia-os/isaura.git

```
- Change the Directory
```
cd isaura

```
- install with pip (use -e for developer mode)
```
pip install -e .

```
*For Windows 10 and above
We recommend [Ubuntu](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) using WSL

Open a Powershell with Admin permissions and type:
```
wsl --install

```
Now, install the [Ubuntu](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) terminal on Windows.
Inside the Ubuntu terminal, you can now follow the installation instructions [above](https://github.com/ersilia-os/ersilia###Install-on-Linux-and-MacOSX)

- Install dependency
```
sudo apt install build-essential

```

- Set up miniconda environment
```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh

```
- Install GitHub CLI
```
conda install gh -c conda-forge

```
- Login with GitHub CLI
```
gh auth login

```
Then follow the instructions
- Install GitHub LFS
```
conda install git-lfs -c conda-forge

```
- Activate GitHub LFS
```
git-lfs install

```
<p align="right">(<a href="#top">back to top</a>)</p>

# Operation
When you finish installing Ersilia, **browse models** in the [Ersilia Model Hub](https://airtable.com/shrXfZ8pqro0jjcsG/tblZGe2a2XeBxrEHP/viwd5XJVLslkE11Tg).

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

You can check the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailed explanations.
<p align="right">(<a href="#top">back to top</a>)</p>

# Contribute
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
- If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
- If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
- If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.
<p align="right">(<a href="#top">back to top</a>)</p>


# Roadmap
We are working to grow the Hub organically and responding to our users' needs. Here is a detail of the next features to come, stay tuned!
[ ] Deployment for Windows System (expected: February 2022)
[ ] Automated third-party model contributions (expected: March 2022)
[ ] Possibility to run lite models online (expected: May 2022)
<p align="right">(<a href="#top">back to top</a>)</p>

# License and Citation
Ersilia is an open-sourced under the GNU Affero General Public License v3.0
<br>Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.
<p align="right">(<a href="#top">back to top</a>)</p>

# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https:/opencollective.com/ersilia).
<p align="right">(<a href="#top">back to top</a>)</p>
