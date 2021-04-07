# Welcome to Ersilia!
[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![EOSI](https://circleci.com/gh/ersilia-os/ersilia.svg?style=svg)](https://circleci.com/gh/ersilia-os/ersilia)
[![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)

This is work in progress.

* Models can be checked at [Ersilia model hub](https://ersilia.io/hub) (coming soon)
* High-level documentation is available at [Ersilia docs](http://ersilia-hub.netlify.app/docs/)
* Low-level documentation is available at [Ersilia Read The Docs](https://ersilia-os.github.io/ersilia/) (work in progress)

## Install

### Mac and Linux

We recommend working inside a [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) environment.
```
conda create -n ersilia python=3.7
conda activate ersilia
```

Then, simply install with pip.
```
pip install ersilia
```

You are done!
```
ersilia --help
```

# Windows Installation

Coming soon...

## Quick start

First, download your model of interest. It will be stored locally at `~/eos`.

```
ersilia fetch eos0abc
```

You can check the catalog of models available in your computer.

```
ersilia catalog --local
```

For more information, check the model card.
```
ersilia card eos0abc
```

Serve your model. A URL will be displayed, together with the APIs available for the model.
```
ersilia serve eos0abc
```

To run the model, use your API of choice. For example, `predict` for Caffeine.
```
ersilia api eos0abc predict "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
```

Don't forget to close the model when you are finished.
```
ersilia close eos0abc
```

If you don't want to use this model anymore, remove it from your computer.
```
ersilia delete eos0abc
```
