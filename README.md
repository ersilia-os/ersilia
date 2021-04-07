# Welcome to Ersilia!
[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![EOSI](https://circleci.com/gh/ersilia-os/ersilia.svg?style=svg)](https://circleci.com/gh/ersilia-os/ersilia)
[![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/)

This is work in progress.

* Models can be checked at [Ersilia model hub](https://ersilia.io/hub) (coming soon)
* High-level documentation is available at [Ersilia docs](http://ersilia-hub.netlify.app/docs/)
* Low-level documentation is available at [Ersilia Read The Docs](https://ersilia-os.github.io/ersilia/) (work in progress)

## Installation

Ersilia can be installed as a [PyPi package](https://pypi.org/project/ersilia/).

```
pip install ersilia
```

To get install the latest version of Ersilia, you can install directly from GitHub.

```
pip install git+git://github.com/ersilia-os/ersilia.git
```

Check CLI.
```
ersilia --help
```

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
