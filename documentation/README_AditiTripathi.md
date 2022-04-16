# Welcome to the Ersilia Model Hub!

![Ersilia_Plum](https://user-images.githubusercontent.com/84270150/162484838-face6352-5508-4156-97ed-dbfcc0ab6249.png)


[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

## Table of Contents:

- [Abouts Us](#about-ersilia)
- [Project Description](#what-is-ersilia-model-hub)
- [Getting started](#installation-guide)
- [Contribution Guidelines](#contribution-guidelines)
- [Future Scope](#future-scope)
- [License and Citation](#license-and-citation)
- [Contact](#contact-us)


# About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io)  is a non-profit organization committed to infectious disease research. Using artificial intelligence techniques, we want to accelerate the discovery of new medications.

# What is Ersilia Model hub?

A free, **open-source** repository of Artificial Intelligence and Machine Learning (**AI/ML**) models for **drug discovery**. Our platform is aimed at helping researchers identify drug candidates for orphan and neglected diseases, design molecules de novo, understand mechanisms of action or anticipate adverse side effects. The ultimate goal of Ersilia is to lower the barrier to drug discovery, encouraging academic groups and companies to pursue the development of new medicines following the principles of **Open Science**. Ersilia disseminates AI/ML models existing in the literature, as well as an in-house collection of models focused on **diseases** that are currently **neglected** by the pharmaceutical industry due to estimated low return on investment. Endemic diseases of low- and middle-income countries (LMIC) belong to this category.

Scientists can use the [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP) to get pre-trained Artificial Intelligence models so they can incorporate these new technologies into their daily experiments.


# Installation Guide

Follow the installation instructions from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/).

Once Ersilia is installed, you can browse models in the Ersilia Model Hub.

## Install on Linux and MacOSX

**Open a terminal. The best is to set up a Conda environment.**
```
 create a conda environment
 
conda create -n ersilia python=3.7
```


```
activate the environment

conda activate ersilia
```



**Then, simply install the Ersilia Python package.**
```
clone from github

git clone https://github.com/ersilia-os/ersilia.git
```
```
cd ersilia
```

```
install with pip (use -e for developer mode)

pip install -e .
```



**You should be done! Quickly check that the CLI works on your terminal.**
```
see ersilia CLI options

ersilia --help
```

## Install on Windows

**We recommend that you install Ubuntu on your Windows machine. This can be now done very easily with WSL. You will need at least Windows 10.**
```
Open a Power Shell with Admin permissions and type:

wsl --install
```
# Contribution Guidelines
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
- If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
- If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
- If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) is followed by the Ersilia Open Source Initiative.


# Future Scope
We are still in the process of expanding the hub to meet the needs of our end users. The following are some of the features we plan to incorporate in the future:

1. Windows-based System Deployment (expected: February 2022)



2. Third-party model contributions that are automated (expected: March 2022)



3. Online running of lite models is possible (expected: May 2022)

# License and Citation
The *MIT License* is used to open-source this repository. If you use it, please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) us.

# Contact us

You can support us via [Open Collective](https:/opencollective.com/ersilia).
- [Twitter](https://twitter.com/ersiliaio) 
- [LinkedIn](https://www.linkedin.com/company/ersiliaio/) 
- [Medium](https://medium.com/ersiliaio)
