<div id="top"></div>
<img src="https://raw.githubusercontent.com/ersilia-os/ersilia/master/assets/Ersilia_Plum.png" height="70">

# 💊 Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Conda Version](https://img.shields.io/conda/vn/conda-forge/ersilia.svg)](https://anaconda.org/conda-forge/ersilia) [![Python 3.8](https://img.shields.io/pypi/pyversions/ersilia
)](https://www.python.org/downloads/release/python-380/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/277068989.svg)](https://zenodo.org/badge/latestdoi/277068989)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)


## Table of Contents

1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Quick Start Guide](https://github.com/ersilia-os/ersilia#quick-start-guide)
3. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
4. [License and Citation](https://github.com/ersilia-os/ersilia#license-and-citation)
5. [About Us](https://github.com/ersilia-os/ersilia#about-us)

## Project Description

The [Ersilia Model Hub](https://ersilia.io) is a unified platform of pre-trained AI/ML models for 🦠 infectious and neglected disease research. Our mission is to offer an open-source, 🛠 low-code solution that provides seamless access to AI/ML models for 💊 drug discovery. Models housed in our hub come from two sources:

- Published models from literature (with due third-party acknowledgement)
- Custom models developed by the Ersilia team or our valued contributors

In Ersilia, you can find models related to antibiotic activity prediction, ADMET prediction, molecular representation, generative chemistry, and much more.

- To find your models of interest, visit the [Ersilia Model Hub browser](https://ersilia.io/model-hub)
- For high-level documentation, check the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
- For a low-level documentation for developers, see the [Ersilia API reference](https://ersilia-os.github.io/ersilia/)

## Quick Start Guide

### Installation

Please check the package requirements in the [Installation Guide](https://ersilia.gitbook.io/ersilia-book/quick-start/installation). In brief:

1. You need to install the Ersilia CLI. This is the central tool to manage models in the Ersilia Model Hub.
2. Each model is packaged independently in its own environment.

To **install the Ersilia CLI**, create a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html) environment and activate it:

```bash
conda create -n ersilia python=3.10
conda activate ersilia
```

Then, clone the current repository and install with `pip`:

```bash
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
pip install -e .
```

Alternatively, you can directly install [from PyPi](https://pypi.org/project/ersilia/):
```bash
pip install ersilia
```

Check that Ersilia is properly installed:
```bash
ersilia --help
```

To **install each model**, you have two options. You can either package them from source, in which case only Conda is required, or as [Docker containers](https://www.docker.com/). We highly recommend having Docker installed and active in your computer. In any case, the Ersilia CLI takes care of the packaging seamlessly.

### Running a Model

Once the Ersilia package is installed, you can use the CLI to run any model seamlessly. First, select a model from the [Ersilia Model Hub](https://ersilia.io/model-hub/) and fetch it. Each model in the Ersilia Model Hub has an associated identifier. For example, the broad spectrum activity prediction model from [_Stokes_ et al. (2020)]() has the identifier `eos4e40`. You can **fetch** it as follows:

```bash
ersilia fetch eos4e40
```

Now you can **serve** the model:

```bash
ersilia serve eos4e40
```

To view some **information** about the model, type the following:

```bash
ersilia info
```

The simplest way to run a model is by passing a CSV file as input. If you don't have one, you can generate it easily. In this case, we take 5 molecules as an **example**:

```bash
ersilia example -n 5 -f my_input.csv
```

Now you can **run** the model:

```bash
ersilia run -i my_input.csv -o my_output.csv
```

To stop the service, you can simply **close** the model:

```bash
ersilia close
```

Finally, if you don't want to use the model anymore, **delete** it as follows:

```bash
ersilia delete eos4e40
```

### Basic Commands

Below is a list of the most important commands of the Ersilia CLI:

| Command    | Description                                      |
|------------|--------------------------------------------------|
| `catalog`  | List a catalog of models                        |
| `close`    | Close model                                     |
| `delete`   | Delete model from local computer                |
| `example`  | Generate input examples for the model of interest |
| `fetch`    | Fetch model from the Ersilia Model Hub          |
| `info`     | Get model information                         |
| `run`      | Run a served model                              |
| `serve`    | Serve model                                     |
| `test`     | Test a model                                    |


Please see the a full reference of all commands available [here](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/developer-docs/command-line-interface).

## Contribute

The Ersilia Model Hub is a free, open source software and we highly value new contributors. There are several ways in which you can contribute to the project:

* A good place to start is checking open [GitHub Issues](https://github.com/ersilia-os/ersilia/issues)
* If you have identified a bug in the code, please open a new issue using the [Bug Report template](https://github.com/ersilia-os/ersilia/issues/new?template=bug_report.yml)
* Share any feedback with the community using [GitHub Discussions](https://github.com/ersilia-os/ersilia/discussions) for the project
* If you wish to contribute to our main codebase, please make sure to follow our [guidelines for codebase quality and consistency](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/developer-docs/developer-guide-for-codebase-quality-and-consistency)
* Check our [Contributing Guide](https://github.com/ersilia-os/ersilia/blob/master/CONTRIBUTING.md) for more details

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) code of conduct.

### Submit a New Model

If you want to incorporate a new model in the platform, open a new issue using the [Model Request template](https://github.com/ersilia-os/ersilia/issues/new?assignees=&labels=new-model&template=model_request.yml&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) or contact us using [this form](https://www.ersilia.io/request-model).

After submitting your model request via an issue (suggested), an Ersilia maintainer will review your request. If they approve your request, a new model respository will be created for you to fork and use following the [eos-template](https://github.com/ersilia-os/eos-template). There is a [demo repository](https://github.com/ersilia-os/eos-demo) explaining the steps one-by-one.

## License and Citation

This repository is open sourced under the [GPL-3 License](https://github.com/ersilia-os/ersilia/blob/master/LICENSE).

If you find this software useful, please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff):

**Ersilia Model Hub: a repository of AI/ML models for infectious and neglected tropical diseases**  
Gemma Turon, Abel Legese, Dhanshree Arora, Miquel Duran-Frigola  
Version: 0.1.43  
DOI: [10.5281/zenodo.7274645](https://doi.org/10.5281/zenodo.7274645)  
GitHub: [https://github.com/ersilia-os/ersilia](https://github.com/ersilia-os/ersilia)

### Authorship

Please note that Ersilia distinguises between software contributors and software authors. The Ersilia Model Hub Authorship guidelines can be found in the [Authorship file](https://github.com/ersilia-os/ersilia/blob/master/AUTHORSHIP.md) and current authors can be found in the [Citation file](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff). We acknowledge past authors of the software below:
- Carolina Caballero

### Cited By

The Ersilia Model Hub is used in a number of scientific projects. Read more about how we are implementing it in:

- [Turon, Hlozek, et al. _Nature Communications_, 2023](https://www.nature.com/articles/s41467-023-41512-2)
- [Van Heerden, et al. _ACS Omega_, 2023](https://pubs.acs.org/doi/10.1021/acsomega.3c05664)
- [Offensperger, Tin, Duran-Frigola, et al. _Science_, 2024](https://www.science.org/doi/10.1126/science.adk5864)
- [Turon, et al. _ACS Med Chem Lett_, 2024](https://doi.org/10.1021/acsmedchemlett.4c00131)

In addition, our views on how the Ersilia Model Hub can be deployed effectively in Africa can be found in these two articles:

- [Turon and Duran-Frigola, _Artificial Intelligence in the Life Sciences_, 2024](https://www.sciencedirect.com/science/article/pii/S2667318524000254?via%3Dihub)
- [Turon and Duran-Frigola, _ACS Infectious Disease Research_, 2024](https://pubs.acs.org/doi/10.1021/acsinfecdis.4c00585)

To see a full list of all Ersilia publications, please visit [this link](https://www.ersilia.io/publications).

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization with the mission to equip laboratories, universities, and clinics in the Global South with AI/ML tools for infectious disease research. We work on the principles of open science, decolonized research, and egalitarian access to knowledge and research outputs. You can support Ersilia by clicking [here](https://www.ersilia.io/donate).

### Funding

The Ersilia Model Hub is our flagship project. The tool is funded via a combination of sources. Full disclosure can be found in our [website](https://ersilia.io/supporters). Highlighted supporters include the [Mozilla Builders Accelerator](https://builders.mozilla.org/), [Fast Forward](https://www.ffwd.org/), [Splunk Pledge](https://www.splunk.com/en_us/about-us/splunk-pledge/nonprofit-license-application.html), the [AI2050 Program by Schmidt Sciences](https://ai2050.schmidtsciences.org/), and the [Spanish Ministry of Science, Innovation, and Universities](https://www.aei.gob.es/convocatorias/buscador-convocatorias/proyectos-generacion-conocimiento-2023) (grant PID2023-148309OA-I00 funded by MICIU/AEI/10.13039/501100011033).

<div id="bottom"></div>
<img src="https://raw.githubusercontent.com/ersilia-os/ersilia/master/assets/ministerio_ciencia_innovacion_universidades_aei.png" height="70">