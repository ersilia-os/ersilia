# Welcome to the Ersilia Open Source Initiative (EOSI) Model Hub! <a name="Top"></a> #

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)

### Table of Contents:
1. [Who We Are](#who-we-are)
2. [Project Description](#description)
3. [Getting started](#getting-started)
4. [Contribute](#contribution)
5. [Roadmap](#roadmap)
6. [License and citation](#license-and-citation)
7. [About us](#about-us)

## Who We Are <a name="who-we-are"></a> ##
The EOSI is an organization that believes that scientific development is a driving component towards the overall progress of Low and Low Middle Income Countries (LMICs).

## Project Description <a name="description"></a> ##
The model hub is our main project and the aim is to provide a user-friendly platform for the deployment of AI/ML models. On this platform, scientists can browse through the assets, identify those which are relevant to their research and obtain predictions without the need of coding expertise. Currently, most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who does not have the necessary expertise. This gap becomes even larger in LMIC institutions where access to bioinformatic facilities or data science experts are scarce. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the licenses of each model before using them.
* Read more about the project in the <a href="https://ersilia.gitbook.io/ersilia-book/">Ersilia gitbook</a>
* Available models can be checked at the <a href = "https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP">Ersilia Model Hub</a>

<div align = "right"><a href="#Top">Back to Top</a></div>

## Getting started<a name="#getting-started"></a> ##
To use Ersilia, kindly follow the installation instructions <a href = "https://ersilia.gitbook.io/ersilia-book/quick-start/installation">here</a>. If this is your first time using Ersilia, once it is installed, let's run through a quick example to familiarize ourselves with how it works.

Browse through the models in the <a href = "https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP">Ersilia Model Hub</a>.

Select one model. For example `chemprop-antibiotic`. You can **fetch** your model with the Ersilia CLI:

```
ersilia fetch chemprop-antibiotic
```
_NB: In case, you decide to go with a different model, the name you use to fetch your model in the Ersilia CLI is the **slug** name on its card._

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

Please see the <a href = "https://ersilia.gitbook.io/ersilia-book/">Ersilia Book</a> for more examples and detailed explanations.

<div align = "right"><a href="#Top">Back to Top</a></div>

## Contribute <a name="#contribution"></a> ##
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:

* If you are a developer, check the <a href = "https://github.com/ersilia-os/ersilia/issues">issues</a> and help us to improve the tool.
* If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please <a href = "https://www.ersilia.io/feedback">contact us</a> or open an issue. We are currently working on an automated contribution tool to facilitate the process.
* If you are a scientist with a cool dataset, also contact us or open an <a href = "https://github.com/ersilia-os/ersilia/issues">issues</a> as we might be interested in developing a model based on the data!
* If there is a third-party model you have identified and would like to see it in the Hub, open an <a href = "https://github.com/ersilia-os/ersilia/issues">issues</a> with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the Contributor Covenant guidelines.

<div align = "right"><a href="#Top">Back to Top</a></div>

## Roadmap <a name="#roadmap"></a> ##
We are working to grow the Hub organically and responding to our users needs. Here are the details of the next features to come, stay tuned!

* Deployment for Windows System (expected: February 2022)
* Automated third-party model contributions (expected: March 2022)
* Possibility to run lite models online (expected: May 2022)
<div align = "right"><a href="#Top">Back to Top</a></div>

## License and citation <a name="#license-and-citation"></a> ##
This repository is open-sourced under the MIT License. Please <a href = "https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff">cite us</a> if you use it.

<div align = "right"><a href="#Top">Back to Top</a></div>

## About us <a name="#about-us"></a> ##
The <a href = "https://ersilia.io/"> Ersilia Open Source Initiative</a> is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via <a href = "https://github.com/opencollective.com/ersilia">Open Collective</a>.

<div align = "right"><a href="#Top">Back to Top</a></div>
