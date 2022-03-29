pp# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

![logo](https://softr-prod.imgix.net/applications/5ce288d5-4600-42df-9d5f-8617b023a3e3/assets/f7b7b7e7-8f53-425e-8504-f722cfa0f503.png)

# Table of Contents:
1. [Project Description](#project-description)
2. [Installation](#installation)
3. [Getting started](#getting-started)
3. [Future development](#future-development)
4. [Contribute](#contribute)
5. [License and citation](#license-and-citation)
6. [About us](#about-us)

# Project Description
[Ersilia Open Source Initiative](https://ersilia.io) is a nonprofit organisation that is focused on research into infectious diseases. We intend to fast-track the development of new medicines by utilizing artificial intelligence tools.

The [Ersilia Model Hub](https://airtable.com/shrUcrUnd7jB9ChZV) is the key project of the Ersilia Open Source Initiative. It is a free repository consisting of Artificial Intelligence and Machine Learning (AI/ML) models. The aim is to provide user friendly models to scientists to aid their research and predict outcomes easily without requiring complex coding expertise.

<p align='justify'>
The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgment) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.
</p>

* Read more about the project at [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* You can check out the available models at the [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)

#
# Installation
To set up the project on your local machine, you need to install a few third-party dependencies before proceeding with the installation of the Ersilia tool. 

They include:-
| <b></b>     | <b></b>           
| :------------------------ | :---------------------------- |
| **`Python`** |
| **`Anaconda or Miniconda (optional)`**            |
| **`Docker (Recommended)`**            |
| **`Git or GitCLI`**            |
| **`Isaura data lake (Recommended)`**      |


The steps for installing this project on a Linux or MacOSX are as follows: 

- First install Python 3.7 upwards from the official [Python website](https://www.python.org/) to your local device

- Install Conda environments, either [anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [miniconda](https://docs.conda.io/en/latest/miniconda.html), for handling dependencies effectively. 
```
- create a conda environment
conda create -n ersilia python=3.7

- activate the environment
conda activate ersilia
```
Detailed instructions for setting up your conda environment can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

- Install [Docker](https://www.docker.com/) to ensure that all the AI/ML models will work smoothly in your local device.

- To enable the Ersilia tool retrieve our codes hosted on [Github](https://github.com/ersilia-os), install [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) or [GitCLI](https://cli.github.com/manual/installation). 

- Clone the ersilia repo.
```
- clone the repo from github
$ git clone https://github.com/ersilia-os/ersilia.git

- change directory into the isausa folder    
$ cd ersilia
      
-  install with pip (use -e for developer mode)
$ pip install -e .
```


- Clone the [Isaura data lake](https://github.com/ersilia-os/isaura) repo. With Isaura, you will be able to store your model predictions on your local computer.

```
- clone the repo from github
$ git clone https://github.com/ersilia-os/isaura.git

- change directory into the isausa folder    
$ cd isaura
      
-  install with pip (use -e for developer mode)
$ pip install -e
```


For more detailed instructions, specific to your operating enviornment, please follow the **installation instructions** from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation).

#
# Getting Started
Once Ersilia is successfully installed, you can search for your AI/ML models of interest in the [Ersilia Model Hub](https://airtable.com/shrXfZ8pqro0jjcsG/tblZGe2a2XeBxrEHP/viwd5XJVLslkE11Tg). All models available in the catalog are accessible through the [Ersilia Python package](https://github.com/ersilia-os/ersilia). 

The image below provides an overview of the Ersilia Model Hub

![systemic-diagram-of-elsilia-model-hub](https://2591732297-files.gitbook.io/~/files/v0/b/gitbook-legacy-files/o/assets%2F-Mj44wxA7bU1hQH19m8I%2F-MjPA_JIOzuKGPtPDw1U%2F-MjPDmAAJkYPLJmKKSD1%2FErsilia_Hub-02.png?alt=media&token=8a876edc-c02e-400c-80c5-b6f3c4060c21)

Select one model. For example `chemprop-antibiotic`. You can **retrieve** your model with the Ersilia CLI:
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

#
# Future development
We are working to grow the Hub organically and responding to our users needs. Here a detail of the next features to come, stay tuned!
1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)

#
# Contribute
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
- If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
- If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
- If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

All **`suggestions`** are welcome!



# License and citation
This repository is open-sourced under the MIT License.
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.

#
# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https:/opencollective.com/ersilia).

For more enquiries, contact us on our [website](https://ersilia.io) or via [email](hello@ersilia.io)
