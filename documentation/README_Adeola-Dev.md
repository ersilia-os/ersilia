# The Ersilia Model Hub

![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)


# Project Description
The Ersilia Model Hub is the main project of the [Ersilia Open Source Initiative](https://ersilia.io). The aim is to provide a platform for a user-friendly deployment of AI/ML models, where scientists can browse through the assets, identify those which is relevant to their research and obtain predictions without the need of coding expertise. Currently, most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who does not have the necessary expertise. This gap becomes even larger in Low and Middle Income Country institutions where access to bioinformatic facilities or data science experts are scarce. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.

* Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)


### Table of content
1. [Getting started](https://github.com/ersilia-os/ersilia#getting-started)
2. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
3. [Roadmap](https://github.com/ersilia-os/ersilia#roadmap)
4. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
5. [About us](https://github.com/ersilia-os/ersilia#about-us)
6. [Our socials](https://github.com/ersilia-os/ersilia#our-socials)



# Getting started

### Installation guide

The following softwares are a prerequisite for the installation of the Erslia Model Hub
1. [Python version 3.7 and above](https://www.python.org/)
2. Conda environment, the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) installers are preferred - they can both be used by the Ersilia project to manage dependencies.
3. [Docker](https://www.docker.com/) - Docker ensures that all our AI/ML assets will work smoothly in your local device.
4. [Github CLI](https://cli.github.com/manual/installation)
5. [Git LFS](https://git-lfs.github.com/)
6. [Ubuntu](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab)  - this is recommended for windows PC.

Follow the **installation instructions** from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation) for full details.

### Testing and running the project models
Once Ersilia is installed, you can **browse models** in the [Ersilia Model Hub](https://airtable.com/shrXfZ8pqro0jjcsG/tblZGe2a2XeBxrEHP/viwd5XJVLslkE11Tg).

One of the models is the **Antibiotic activity prediction Model**. In this case example, we show how to run predictions based on the AI/ML model.

The example `chemprop-antibiotic` will be used in this case. You can **fetch** your model with the Ersilia CLI - this command installs and downloads the model to your local PC
```
ersilia fetch chemprop-antibiotic
```
Generate a few (5) example molecules, to be used as input. Molecules are typically expressed in SMILES format.
```
ersilia example chemprop-antibiotic -n 5 -f my_molecules.csv
```
Once the model has been fetched, it should be ready to be used. A model in the Ersilia Model Hub can be thought of as a set of APIs. You can serve the model like this:
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

# Contribute
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
* If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
* If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
* If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
* If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

A special credit goes to;

[Gemma Turon](https://github.com/GemmaTuron)  and 

[Miquel Duranfrigola](https://github.com/miquelduranfrigola) for their work.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

# Roadmap
We are working to grow the Hub organically and responding to our users needs. Here are the details of the next features to come, stay tuned!

1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)

# License and citation
This repository is open-sourced under the GNU AFFERO GENERAL PUBLIC LICENSE. Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.

# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https:/opencollective.com/ersilia).

# Our Socials
Join the Ersilia contributor circle on [Slack Workspace](https://ersilia-workspace.slack.com/)

Follow Ersilia on [Twitter](https://twitter.com/ersiliaio)
