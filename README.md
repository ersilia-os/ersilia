# Welcome to Ersilia!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia-hub.netlify.app/docs/) [![EOSI](https://circleci.com/gh/ersilia-os/ersilia.svg?style=svg)](https://circleci.com/gh/ersilia-os/ersilia) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Github All Releases](https://img.shields.io/github/downloads/ersilia-os/ersilia/total.svg)](./) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

* Models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP) \(in preparation\)
* High-level documentation is available at [Ersilia docs](http://ersilia-hub.netlify.app/docs/)
* Low-level documentation is available at [Ersilia Read The Docs](https://ersilia-os.github.io/ersilia/)
* Learn about us at: [Ersilia Website](https://ersilia.softr.app)

## Install

### Mac and Linux

We recommend working inside a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) environment.

```text
conda create -n ersilia python=3.7
conda activate ersilia
```

Then, simply install with pip.

```text
pip install ersilia
```

You are done!

```text
ersilia --help
```

### Windows Installation

Coming soon...

## Quick start

First, download your model of interest. It will be stored locally at `~/eos`.

```text
ersilia fetch chemprop-antibiotic
```

You can check the catalog of models available in your computer.

```text
ersilia catalog
```

For more information, check the model card.

```text
ersilia card chemprop-antibiotic
```

Serve your model. A URL will be displayed, together with the APIs available for the model.

```text
ersilia serve chemprop-antibiotic
```

To run the model, use your API of choice. For example, `predict` for Caffeine.

```text
ersilia api chemprop-antibiotic predict "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
```

Don't forget to close the model when you are finished.

```text
ersilia close chemprop-antibiotic
```

If you don't want to use this model anymore, remove it from your computer.

```text
ersilia delete chemprop-antibiotic
```

# Special thanks

We would like to thank the Open source community and, especially, [FossHost](https://fosshost.org/) for providing help with infrastructure.
