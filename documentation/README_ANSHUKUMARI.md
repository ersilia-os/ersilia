# Welcome to Ersillia Model Hub!







[![PayPal](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012)
[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![AGPL License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0)
[![Contributor](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](https://github.com/ersilia-os/ersilia/blob/master/documentation/code_of_conduct.md)
[![Passed](https://circleci.com/gh/ersilia-os/ersilia.svg?style=svg)](https://circleci.com/gh/ersilia-os/ersilia)
[![pypipackage](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/)
[![python](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![Gitpod](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)
[![Maintained](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/Naereen/StrapDown.js/graphs/commit-activity)
[![linux](https://svgshare.com/i/Zhy.svg)]()
[![macOS](https://svgshare.com/i/ZjP.svg)]()
[![](https://svgshare.com/i/ZhY.svg)]()




## Table Of Content

1.	Project Description
2.  Overview
3.	Technologies
4.	Contribution Guidelines 
5.	Prerequisites
6.	Features
7.	Installation Instructions
8. RoadMap
9.	Requirements
10.	Operation 
11.	License
12.	Vision
13.	Future Scope 
14.	Attributes
15.	Acknowledgements
16.	Credits/Contacts
17.	About Us

## Project Description

The Ersilia Model Hub is a platform which provide you a pre trained Artificial Intelligence model for scientists so that the can integrate it in daily tasks. Ersilia Open Source Initiative is a part of the Ersilia Model Hub. It is creating a software known as FOSS software and this software should be accompanied by a good documentation, usage, contribution, instruction of how to install and to use and many more about it. FOSS software which basically focuses to create a platform on AI and ML models for finding assets relevant to their research and can try them without coding expertise. 

This project is inspired by the participation of the Ersilia Open Source Initiative in three programs this year: The [Digital Infrastructure Incubator](https://github.com/matiassingers/awesome-readme), focused on FOSS community and governance, the [Open Life Sciences](https://openlifesci.org/), , for implementing Open Science practices at all levels and the [Software Sustainability Fellowship](https://www.software.ac.uk/programmes-and-events/fellowship-programme), oriented towards the sustainability of research software tools. The intern will gain access to all documentation and lessons learned during these programs and help implement them.

## Overview

![.](https://user-images.githubusercontent.com/63330165/160850205-9d269457-06ad-46b7-9aaa-7a934c2fb47c.png)


## Technologies
## Contribution Guidelines

The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:

1. If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
2. If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://www.ersilia.io/) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
3. If you are a scientist with a cool dataset, also [contact us](https://www.ersilia.io/) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
4. If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

## Prerequisites

Please make sure you have the installation of below libraries on your machine.
1.	[Python 3.7 and above](https://www.python.org/)
2.	[Docker](https://www.docker.com/)
3.	[Conda](https://www.docker.com/) / [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) / [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
4.	[Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [Git LFS](https://git-lfs.github.com/)

## Features

| Feature           | Description                                                              |
| ----------------- | ------------------------------------------------------------------ |
| Bioactivity signatures | Tools are designed to facilitate the use of AI/ML tools. Anyone even from non-coding background can browse a collection of models, choose the ones relevant to their research interests and run predictions without writing even a single line of code. |
| Data driven | Ersilia technology achieves state-of-the-art performance thanks to the integration of chemical, genomic and biomedical text data.|
| User friendly platform | Our tools are designed to facilitate the use of AI/ML tools and platform is user friendly, it is beginners-friendly as well. |
| Open source | All our assets are released under the MIT license. This means the anyone in ersilia community can review, contribute and improve our code, resulting in tools validated more extensively than in the traditional peer-review system. |


##  Installation Instruction

### For Windows
We recommend that you install Ubuntu on your Windows. This can be now done very easily with WSL. You will need at least Windows 10.

	Open a Power Shell with Admin permissions and type:
```bash
  1.	wsl –install
```


	Then simply install the Ubuntu terminal on Windows.

	Inside the Ubuntu terminal. Open a terminal. The best is to set up a Conda environment.
```bash
1.    # create a conda environment
2.	  conda create -n ersilia python=3.7
3.	  # activate the environment
4.    conda activate ersilia
```


	Then, simply install the Ersilia Python package.
```bash
1.	  # clone from github
2.    git clone https://github.com/ersilia-os/ersilia.git
3.	  cd ersilia
4.    # install with pip (use -e for developer mode)
5.    pip install -e .
```


	You should be done! Quickly check that the CLI works on your terminal.
```bash
1.    # see ersilia CLI options
2.	  ersilia --help

```


### For MacOSX and Linux

	Open a terminal. The best is to set up a Conda environment.
```bash
1.	  # create a conda environment
2.	  conda create -n ersilia python=3.7
3.	  # activate the environment
4.	  conda activate ersilia

```


	Then, simply install the Ersilia Python package.
```bash
1.	  # clone from github
2.    git clone https://github.com/ersilia-os/ersilia.git
3.    cd ersilia
4.    # install with pip (use -e for developer mode)
5.    pip install -e .

```


	You should be done! Quickly check that the CLI works on your terminal.
```bash
1.    # see ersilia CLI options
2.	  ersilia --help

```

## RoadMap

We are working to grow the Hub organically and responding to our users needs.
     Here a detail of the next features to come.
1.	Deployment for Windows System (expected: February 2022)
2.	Automated third-party model contributions (expected: March 2022)
3.	Possible to run lite models online (expected: May 2022)

For more clarification visit: https://ersilia.gitbook.io/ersilia-book/quick-start/installation

## Requirements
Note – Do check the licenses of every model before use 

bentoml==0.11.0

PyGitHub

streamlit

pygit2

osfclient

joblib

hashids

bioservices

biopython

markdown

pycrypto

python-dateutil<2.8.1,>=2.1

urllib3

requests<=2.24

dockerfile-parse

pytest==3.10.0

sphinx

boto3

PyYAML

emoji

loguru

virtualenv

pyairtable

PyDrive2

h5py

inputimeout

## Operations

The Ersilia Open Source Initiative is a UK-based charity focused on strengthening the research capacity for infectious and neglected diseases by developing and implementing novel artificial intelligence and machine learning tools.
ML is advancing at incredible speed, holding the promise to revolutionize biomedical practice by leveraging large amounts of data and decreasing the number and cost of laboratory experiments. Unfortunately, in practice, many of ML benefits are only accessible to experienced data scientists, without real adoption of the technology in day-to-day research in laboratories and hospitals.

## License

This repository is open-sourced under [MIT License](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff). Please check the link bellow before using it.

## Vision

The project is focused on strengthening the community guidelines, working on materials for and dissemination like blogposts, scientific publications etc. The mission is to again maintain the balance between countries by supporting research in unfunded settings and remove the gap created.

## Future Scope

It will help you to create an effective documentation. The tool developed worldwide are still untouched who doesn’t have the expertise even after they are fully published on open source. Due to these untouched tools a huge gap has been created specially in low and middle income countries. Through this access to bioinformatic facilities or data science experts are scarce.  So to heal this problem we have created FOSS software to ease the access to expertise and support research into neglected infectious diseases.
This Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public spaces. Examples of representing our community include using an official e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event.
## Attributes

This Code of Conduct is adapted from the Contributor Covenant, version 2.0, available at https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.
Community Impact Guidelines were inspired by Mozilla's code of conduct enforcement ladder.
For answers to common questions about this code of conduct, see the FAQ at https://www.contributor-covenant.org/faq. Translations are available at https://www.contributor-covenant.org/translations.
## Acknowledgements

I would like to express my special thanks of gratitude to whole team of outreachy as well as the team of Ersilia and Gemma Turon who gave me this golden opportunity to do this wonderful project based on documentary on behalf of outreachy, which also helped me in learning lot more new things.
Secondly, I would also like a thank my whole team and friends who helped me a lot in creating it.
It helped me increase my knowledge and skills.


## Credit

The credit of this goes to the [Digital Infrastructure Incubator](https://github.com/matiassingers/awesome-readme), focused on FOSS community and governance, the [Open Life Sciences](https://openlifesci.org/), for implementing Open Science practices at all levels and the [Software Sustainability Fellowship](https://www.software.ac.uk/programmes-and-events/fellowship-programme)

## About Us

The [Ersilia Open Source Initiative](https://www.ersilia.io/) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.
You can support us via Open Collective.