<div id="top" align="center">

# Welcome to the **Ersilia Model Hub**!
<br></br>
![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)
<br></br>
[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

</div>
<br></br>

## Table of Contetns:

- [Welcome to the **Ersilia Model Hub**!](#welcome-to-the-ersilia-model-hub)
  - [Table of Contetns:](#table-of-contetns)
  - [Abouts Us](#abouts-us)
  - [Project Description](#project-description)
  - [Features](#features)
  - [Technology](#technology)
  - [Prerequisites](#prerequisites)
  - [Getting started](#getting-started)
  - [Contribution Guidelines](#contribution-guidelines)
  - [Roadmap](#roadmap)
  - [License and citation](#license-and-citation)
  - [Credit](#credit)
  - [Contact](#contact)

<br></br>

## Abouts Us

The [Ersilia Open Source Initiative](https://www.ersilia.io/) is a Non Profit Organization registered with the Charity Commission for England and Wales (no. 1192266). It is our mission to reduce the gap between the productivity of biomedical research between countries by promoting research in underfunded settings and making scientific advances accessible to everyone. You can support us through [Open Collective](https://opencollective.com/ersilia#category-CONTRIBUTE).

<p align="right">(<a href="#top">back to top</a>)</p>

## Project Description

The Ersilia Model Hub is the main project of the [Ersilia Open Source Initiative](https://www.ersilia.io/). It is an open-source repository of Artificial Intelligence and Machine Learning (AI/ML) models for drug discovery. The goal is to provide platforms for the deployment of AI/ML models that scientists can browse through, identify relevant assets and generate predictions without requiring coding expertise. A large majority of the scientific community does not possess the required expertise to use most of the tools developed, even when they are fully open-sourced and published. The gap is even more pronounced in institutions in Low and Middle Income Countries where bioinformatics facilities and data science experts are scarce. Through this project, we aim to democratize access to this expertise and promote research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.

* Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)
<p align="right">(<a href="#top">back to top</a>)</p>

## Features

| Feature           | Description                                                              |
| ----------------- | ------------------------------------------------------------------ |
| Bioactivity signatures | Tools are designed to facilitate the use of AI/ML tools. Anyone even from non-coding background can browse a collection of models, choose the ones relevant to their research interests and run predictions without writing even a single line of code. |
| Data driven | Ersilia technology achieves state-of-the-art performance thanks to the integration of chemical, genomic and biomedical text data.|
| User friendly platform | Our tools are designed to facilitate the use of AI/ML tools and platform is user friendly, it is beginners-friendly as well. |
| Open source | All our assets are released under the GPLv3 license. This means the anyone in ersilia community can review, contribute and improve our code and can distribute their own versions. |

<p align="right">(<a href="#top">back to top</a>)</p>

## Technology

The Ersilia Hub Model uses the Chemical Checker Ersilia's backbone technology developed in patrick aloy lab, IRB Barcelonia. The Chemical Checker encodes chemical information and biological information of the model, it encodes chemical and biological information of this model molecules and represent them in numerical vector so that a computer can analyse.

<p align="right">(<a href="#top">back to top</a>)</p>

## Prerequisites

Please make sure you have the installation of below libraries on your machine.
1.	[Python 3.7 and above](https://www.python.org/)
2.	[Docker](https://www.docker.com/)
3.	[Conda](https://www.docker.com/) / [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) / [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
4.	[Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) and [Git LFS](https://git-lfs.github.com/)

<p align="right">(<a href="#top">back to top</a>)</p>

## Getting started
Follow the **installation instructions** from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation).

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
<p align="right">(<a href="#top">back to top</a>)</p>

## Contribution Guidelines

The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:

1. If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
2. If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://www.ersilia.io/) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
3. If you are a scientist with a cool dataset, also [contact us](https://www.ersilia.io/) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
4. If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

<p align="right">(<a href="#top">back to top</a>)</p>

## Roadmap
We are working to grow the Hub organically and responding to our users needs. Here are the details of the next features to come, stay tuned!
1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)
   
<p align="right">(<a href="#top">back to top</a>)</p>

## License and citation

This repository is open-sourced under GNU General Public License v3.0 or later. See [LICENSE](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) to see the full text. Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff), if you use it.

<p align="right">(<a href="#top">back to top</a>)</p>

## Credit

This project is inspired by the participation of the Ersilia Open Source Initiative in three programs this year: The [Digital Infrastructure Incubator](https://github.com/matiassingers/awesome-readme), focused on FOSS community and governance, the [Open Life Sciences](https://openlifesci.org/), for implementing Open Science practices at all levels and the [Software Sustainability Fellowship](https://www.software.ac.uk/programmes-and-events/fellowship-programme), oriented towards the sustainability of research software tools.

## Contact

You can support us via [Open Collective](https:/opencollective.com/ersilia).
- [Twitter](https://twitter.com/ersiliaio) - [LinkedIn](https://www.linkedin.com/company/ersiliaio/) - [Medium](https://medium.com/ersiliaio)
<p align="right">(<a href="#top">back to top</a>)</p>
