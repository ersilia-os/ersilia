# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) 

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![EOSI](https://circleci.com/gh/ersilia-os/ersilia.svg?style=svg)](https://circleci.com/gh/ersilia-os/ersilia) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Github All Releases](https://img.shields.io/github/downloads/ersilia-os/ersilia/total.svg)](./) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

![logo](https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png)


### Table of Contents:
1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Getting started](https://github.com/ersilia-os/ersilia#getting-started)
3. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
4. [Roadmap](https://github.com/ersilia-os/ersilia#roadmap)
5. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
6. [About us](https://github.com/ersilia-os/ersilia#about-us)

# Project Description
With just some compound inputs and pretrained AI tools, scientists can now have access to millions of data for research. The Ersilia Model Hub, a [project](https://ersilia.gitbook.io/ersilia-book/) of the Ersilia Open Source Initiative made this possible. It embeds the experimental data for each compound, plus targets and side-effects profiles, thus making it a more powerful and clinically relevant predictor. You can see the available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)


## Key Features
* GitHub Flavored Markdown
* Models for infectious and neglected diseases 
* Large scale datasets
* Scientific literature
* Full-screen mode
* Work distraction-free
* Windows and Linux ready. 
 

# Getting started
You should follow the **installation instructions** from the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation).

Once Ersilia is installed, you can **browse models** in the [Ersilia Model Hub](https://airtable.com/shrXfZ8pqro0jjcsG/tblZGe2a2XeBxrEHP/viwd5XJVLslkE11Tg).

Select one model. For example `chemprop-antibiotic`. You can **fetch** your model with the Ersilia CLI:
```
ersilia fetch chemprop-antibiotic
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

Please see the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailes explanations.

# How to Contribute
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! 

If you are a developer who wants to help us to improve the tool, have developed a model, a scientist with a cool dataset, or you want us to incorporate a third-party model in the Hub, open an issue with the relevant information and we will get back to you as soon as possible. Here’s how:

1. Clone this repo https://github.com/ersilia-os/ersilia and create a new branch
2. Make the changes you want 
3. Submit Pull Request with a comprehensive description of the change

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

# Roadmap
We are working to grow the Hub organically and responding to our users' needs. Here a detail of the next features to come, stay tuned!
Outreachy applicants, these are ideas we can work on during your internship.
1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)

# License and citation
This repository is open-sourced under the GNU Affero General Public License.
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.

# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https:/opencollective.com/ersilia).

# Links
* [Website](https://ersilia.io)
* [Documentation](https://github.com/ersilia-os/ersilia/tree/master/documentation)
* [Issue tracker](https://github.com/ersilia-os/ersilia/issues)

