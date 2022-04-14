pp# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

![logo](https://softr-prod.imgix.net/applications/5ce288d5-4600-42df-9d5f-8617b023a3e3/assets/f7b7b7e7-8f53-425e-8504-f722cfa0f503.png)

> # Table of Contents:
1. [Project Description](#project-description)
2. [Elrisia Model Hub](#ersilia-model-hub)
3. [Installation](#installation)
    - [Installation on MacOSX or Linux](#installation-on-macosx-or-linux)
    - [Installation on Windows](#installation-on-windows)
4. [Retrieve Ersilia codes from GitHub](#retrieve-ersilia-codes-from-github)
    - [Clone the Ersilia repo](#clone-the-ersiliahttpsgithubcomersilia-osersiliagit-repo)
    - [Clone the Isaura data lake repo](#clone-the-isaura-data-lakehttpsgithubcomersilia-osisauragit-repo)
5. [Getting started](#getting-started)
    - [Usage Steps](#usage-steps)
6. [Future development](#future-development)
7. [Contribute](#contribute)
8. [License and citation](#license-and-citation)
9. [About us](#about-us)
    

> # Project Description
[Ersilia Open Source Initiative](https://ersilia.io) is a nonprofit organisation that is focused on the research into infectious diseases. We intend to fast-track the development of new medicines by utilizing artificial intelligence tools.

The [Ersilia Model Hub](https://airtable.com/shrUcrUnd7jB9ChZV) is the key project of the Ersilia Open Source Initiative. It is a free repository which contains Artificial Intelligence and Machine Learning (AI/ML) models. The aim is to provide user friendly models to scientists to aid their research and predict outcomes easily without requiring complex coding abilities.

<p align='justify'>
Currently, majority of the developed tools, including published and open-sourced, remain unusable to a vast majority of the scientific community who lack the necessary expertise. This is problematic especially in Low and Middle Income Country institutions with limited access to bioinformatic facilities or data science experts. This project aims to provide access to this expertise and support research into a variety of neglected and infectious diseases.

<p align='justify'>
The models embedded in the hub include established models as well as propriety models developed by the Ersilia team.
</p>

* Read more about the project at [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* You can check out the available models at the [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)

#
> # Ersilia Model Hub

![systemic-diagram-of-elsilia-model-hub](https://2591732297-files.gitbook.io/~/files/v0/b/gitbook-legacy-files/o/assets%2F-Mj44wxA7bU1hQH19m8I%2F-MjPA_JIOzuKGPtPDw1U%2F-MjPDmAAJkYPLJmKKSD1%2FErsilia_Hub-02.png?alt=media&token=8a876edc-c02e-400c-80c5-b6f3c4060c21)

| <b>Attributes</b>     | <b>Description</b>           
| :------------------------ | :---------------------------- |
| **`Open Source`** |All the Ersilia assets are released under the MIT license. This enables the scientific community to review, contribute and improve our code. The advantage is the development of tools which have been validated more extensively than in the usual peer-review system. |
| **`User friendly`**|Our tools are designed with scientists without complex coding expertise in mind, thus they are beginner friendly. Anyone can browse and make a choice from the collection of models based on their research interest and run predictions with ease.  |
| **`Bioactivity Signatures`**|The models available in the hub are built upon a rich representation of molecules consisting of extensive experimental data including targets and side effect profiles, resulting in more clinically relevant and powerful predictions. |
| **`Data driven`**|Ersilia technology achieves a state-of-the-art performance due to the integration of chemical, genomic and biomedical text data. Our AI tools are extensively trained using millions of data from scientific literature.|



#
> # Installation
To set up the project on your local machine, you need to install a few third-party dependencies before proceeding with the installation of the Ersilia tool. 

They include:
| <b></b>     | <b></b>           
| :------------------------ | :---------------------------- |
| [Python](https://www.python.org/) |
| [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) `(Optional)`       |
| [Docker](https://www.docker.com/) `(Recommended)`            |
| [Git](https://git-scm.com/download/win) or [GitCLI](https://github.com/cli/cli#installation)     `(Recommended)`      |
| [Isaura data lake](https://github.com/ersilia-os/isaura.git) `(Recommended)`     |

#
> ## Installation on MacOSX or Linux
The steps for installing this project on a Linux or MacOSX are as follows: 

- First install Python 3.7 upwards from the official [Python website](https://www.python.org/) to your local device

- Install Conda environments, either [anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html), for handling dependencies effectively. 
```
# to create a conda environment
$ conda create -n ersilia python=3.7

# to activate the environment
$ conda activate ersilia
```
Detailed instructions for setting up the conda environment can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

- Install [Docker](https://www.docker.com/). This  ensures that all the AI/ML models will work seamlessly on your local device.

#
> ## Installation on Windows
The Elrisia tool was optimized for a system with Ubuntu OS. For installation on a windows system, we recommend that Ubuntu is installed. This can be done very easily using the Windows Subsystem for Linux (WSL).


WSL is a free, optional feature of Windows 10 that allows Linux programs to run on Windows. 

- To install WSL, run Command Prompt as Administrator, then type the following command.
```
$ wsl --install -d Ubuntu
```
- Once the Ubuntu terminal is up and running, you'll be prompted to create a default UNIX user account. Input a username and password.

- Next download anaconda or miniconda. 

`Note*`

The installation steps outlined below were performed with Anaconda v5.3.1

```
wget https://repo.continuum.io/archive/Anaconda3-5.3.1-Linux-x86_64.sh
```
- Run the Anaconda installation script

```
bash Anaconda3-5.3.1-Linux-x86_64.sh
```
- Follow the prompts and accept the license agreement.

`Note*` 

Confirm that Anaconda and Python were installed by running `'which python'` in the terminal, you should see a path to your installation i.e. `/home/pauline/anaconda3/bin/python`.

Once Anaconda is successfully installed, proceed with the steps to activate the ersilia conda environment.

```
$ conda create -n ersilia python=3.7

# to activate the environment
$ conda activate ersilia
```
#


For more detailed instructions, specific to your operating enviornment, please follow the instructions outlined in the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation).
#

> # Retrieve Ersilia codes from GitHub.
- To enable the Ersilia tool retrieve the codes hosted on [Github](https://github.com/ersilia-os), install [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) or [GitCLI](https://cli.github.com/manual/installation). 


> ## Clone the [Ersilia](https://github.com/ersilia-os/ersilia.git) repo.

This enables you to create a local copy of the Ersilia repo on your computer.
- Clone the repo from github
```
$ git clone https://github.com/ersilia-os/ersilia.git
```
- Change directory into the Ersilia folder 
```
$ cd ersilia
```      
- Install dependencies
```
# Use -e . to install from setup.py using developer mode
$ pip install -e .

OR

# install from .txt file
$ pip install -r requirements.txt
```
#
> ## Clone the [Isaura data lake](https://github.com/ersilia-os/isaura.git) repo.

This enables you to backup your model predictions on your local computer.
- Clone the repo from github
```
$ git clone https://github.com/ersilia-os/isaura.git
```
- Change directory to the isaura folder    
```
$ cd isaura
```      
- Install dependencies
```
# Use -e . to install from setup.py using developer mode
$ pip install -e .

OR

# install from .txt file
$ pip install -r requirements.txt
```


#
> # Getting Started
Once Ersilia is successfully installed, you can search for your AI/ML models of interest in the [Ersilia Model Hub](https://airtable.com/shrXfZ8pqro0jjcsG/tblZGe2a2XeBxrEHP/viwd5XJVLslkE11Tg). All available models are accessible through the [Ersilia Python package](https://github.com/ersilia-os/ersilia). 

`Note*`

You need to use your linux terminal with the activated ersilia environment to run these commands.

> ## Usage Steps
- First fetch the model of interest:
```
$ ersilia fetch <name-of-model>
```
- Generate some example molecules, to be used as input:
```
$ ersilia example <name-of-model> -n 5 -f my_molecules.csv
```
- Serve your model:
```
$ ersilia serve <name-of-model>
```
- Run the prediction API:
```
$ ersilia api -i my_molecules.csv -o my_predictions.csv
```
- Finally, close the service when you are done.
```
$ ersilia close
```
More detailed instructions and examples can be found in the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/).

#
># Future development
We are still undergoing processes in order to expand the hub, taking into account the needs from our end users. Some of the features we intend to implement in the future include:
- Deployment for Windows based Systems (expected: February 2022)

- Automated third-party model contributions (expected: March 2022)

- Possibility to run lite models online (expected: May 2022)

#
> # Contribute

<p aligh='justify'>
The Ersilia Model Hub is constantly being developed and maintained by a team of Ersilia employees and volunteers. If you would like to contribute to this project, there are several ways in which you can!
</p>

- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and assist us in making imporvements to the ersilia tool.

- If you have some experience developing models and would like to include your models in the Hub, [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues).

- If you are a scientist with an interesting dataset, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data.

- If there is a third-party model you have identified and would like us to add it to the hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

All **`suggestions`** are welcome!



> # License and citation
This repository is open-sourced under the GNU Affero General Public License. Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.

#
> # About us

The [Ersilia Open Source Initiative](https://ersilia.io) is incorporated with the Charity Commission for England and Wales (number 1192266). We are passionate about reducing the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https:/opencollective.com/ersilia).

For more enquiries, contact us on our [website](https://ersilia.io) or via [email](hello@ersilia.io).