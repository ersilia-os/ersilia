<div id="top"></div>
<img src="https://raw.githubusercontent.com/ersilia-os/ersilia/master/assets/Ersilia_Plum.png" height="70">

# 🎉 Welcome to the Ersilia Model Hub 🌟

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Conda Version](https://img.shields.io/conda/vn/conda-forge/ersilia.svg)](https://anaconda.org/conda-forge/ersilia) [![Python 3.8](https://img.shields.io/pypi/pyversions/ersilia
)](https://www.python.org/downloads/release/python-380/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black)
[![DOI](https://zenodo.org/badge/277068989.svg)](https://zenodo.org/badge/latestdoi/277068989) [![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)


## Table of Contents

1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Quick start guide](https://github.com/ersilia-os/ersilia#quick-start-guide)
3. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
4. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
5. [About us](https://github.com/ersilia-os/ersilia#about-us)

## Project Description

The [Ersilia Model Hub](https://ersilia.io) is a unified platform of pre-trained AI/ML models for 🦠 infectious and neglected disease research. Our mission is to offer an open-source, 🛠 low-code solution that provides seamless access to AI/ML models for 💊 drug discovery. Models housed in our hub come from two sources:

- Published models from literature (with due third-party acknowledgement)
- Custom models developed by the Ersilia team or our valued contributors

You can read more about the project in the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) and browse available models in the [Ersilia Model Hub](https://ersilia.io/model-hub/).

## Quick Start Guide

Please check the package requirements in the [Installation Guide](https://ersilia.gitbook.io/ersilia-book/quick-start/installation). The following steps are a quick start guide to using Ersilia.

First, create a conda environment and activate it:

```bash
conda create -n ersilia python=3.10
conda activate ersilia
```

Then, clone this repository and install with `pip`:

```bash
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
pip install -e .
```

Alternatively, you can directly install from PyPi:
```bash
pip install ersilia
```

Once the Ersilia package is installed, you can use the CLI to run predictions. First, select a model from the [Ersilia Model Hub](https://ersilia.io/model-hub/) and fetch it:

```bash
ersilia fetch eos4e40
```

Note that you can use the model identifier (eos4e40) or its human-readable slug (antibiotic-activity).

Now you can serve the model:

```bash
ersilia serve eos4e40
```

To view some information of the model, type the following:

```bash
ersilia info
```

The simplest way to run a model is by passing a CSV file as input. If you don't have one, you can generate it easily. In this case, we take 5 molecules as an example:

```bash
ersilia example -n 5 -f my_input.csv
```

Now you can run the model:

```bash
ersilia run -i my_input.csv -o my_output.csv
```

To stop the service, you can simply close the model:

```bash
ersilia close
```

Finally, if you don't want to use the model anymore, delete it as follows:

```bash
ersilia delete eos4e40
```

Please see the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailed explanations.

For Python versions 3.12, Ersilia explicitly installs the setuptools library during installation. This is due to a compatibility issue in Python 3.12, which is described in python/cpython#95299.

Note: If you are using Python 3.12, you don’t need to take any manual action. The Ersilia CLI automatically handles this by installing setuptools as part of the setup process.

## Contribute

The Ersilia Model Hub is a Free, Open Source Software and we highly value new contributors. There are several ways in which you can contribute to the project:

* A good place to start is checking open [issues](https://github.com/ersilia-os/ersilia/issues)
* If you have identified a bug in the code, please open a new issue using the bug template
* Share any feedback with the community using [GitHub Discussions](https://github.com/ersilia-os/ersilia/discussions) for the project
* Check our [Contributing Guide](https://github.com/ersilia-os/ersilia/blob/master/CONTRIBUTING.md) for more details

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) code of conduct.

### Development Guidelines

To maintain consistency and code quality, we follow certain coding and linting standards. Please adhere to these guidelines when contributing:

#### Pre-commit Hooks

We use `pre-commit` and `ruff` to automate code quality checks. Ensure you install and set up `pre-commit` and `ruff` before committing any changes:

1. Install pre-commit: `pip install pre-commit`
2. Set up pre-commit hooks in your local repository by running:
   ```bash
   pre-commit install
   ```
3. When you commit it automatically fix the issues but will fail for critical error such as missing docstring on a public class and public methods.

#### Manual with Ruff

1. Run `ruff` to check for linting errors:
   ```bash
   ruff check .
   ```
2. Automatically fix linting issues (where possible):
   ```bash
   ruff check . --fix
   ```

#### Docstring Style

We adhere to the [NumPy-style docstring format](https://numpydoc.readthedocs.io/en/latest/format.html). Please document all public methods and functions using this style.

Consistent documentation ensures the code is easy to understand and maintain.

Thank you for your contributions and for helping make the Ersilia Model Hub a better project!

### Submit a New Model

If you want to incorporate a new model in the platform, open a new issue using the [model request template](https://github.com/ersilia-os/ersilia/issues/new?assignees=&labels=new-model&template=model_request.yml&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) or contact us using the following [form](https://www.ersilia.io/request-model).

After submitting your model request via an issue (suggested), an Ersilia maintainer will review your request. If they approve your request, a new model respository will be created for you to fork and use! There is a [demo repository](https://github.com/ersilia-os/eos-demo) explaining the steps one-by-one.

## License and Citation

This repository is open-sourced under the [GPL-3 License](https://github.com/ersilia-os/ersilia/blob/master/LICENSE).
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it!

### Authorship

Please note that Ersilia distinguises between software contributors and software authors. The Ersilia Model Hub Authorship guidelines can be found in the [Authorship file](https://github.com/ersilia-os/ersilia/blob/master/AUTHORSHIP.md) and current authors can be found in the Citation file. We acknowledge past authors of the software below:
- Carolina Caballero

### Cited by

The Ersilia Model Hub is used in a number of scientific projects. Read more about how we are implementing it in:
- [Turon, Hlozek et al, Nat Commun, 2023](https://www.nature.com/articles/s41467-023-41512-2)
- [Van Heerden et al, ACS Omega, 2023](https://pubs.acs.org/doi/10.1021/acsomega.3c05664)
- [Offensperger et al, Science, 2024](https://www.science.org/doi/10.1126/science.adk5864)
- [Turon et al, ACS Med Chem Lett, 2024](https://doi.org/10.1021/acsmedchemlett.4c00131)

## About Us

The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.
[Help us](https://www.ersilia.io/donate) achieve our mission!

### Funding

The Ersilia Model Hub is the flagship product of Ersilia. It has been funded thanks to a combination of funding sources. Full disclosure can be found in our [website](https://ersilia.io/supporters). Highlighted supporters include the Mozilla Builders Accelerator, Fast Forward, Splunk Pledge and the AI2050 Program by Schmidt Sciences. 
