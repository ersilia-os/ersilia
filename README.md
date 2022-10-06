<div id="top"></div>
<img src="https://github.com/ersilia-os/ersilia/blob/master/assets/Ersilia_Plum.png" height="70">

# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![License: AGPL v3](https://img.shields.io/badge/License-AGPL%20v3-yellow.svg)](https://www.gnu.org/licenses/agpl-3.0)

[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)

### Table of Contents:
1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Installation](https://github.com/ersilia-os/ersilia#installation)
3. [Contribute](https://github.com/ersilia-os/ersilia#contribute)
4. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
5. [About us](https://github.com/ersilia-os/ersilia#about-us)

# Project Description
The Ersilia Model Hub is a unified platform of pre-trained AI/ML models for infectious and neglected disease research. The end goal is to provide an open-source, no-code solution to access AI/ML models to accelerate drug discovery. The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors.

* Read more about the project in the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Browse available models in the [Ersilia Model Hub](https://ersilia.io/model-hub/)
<p align="right">(<a href="#top">back to top</a>)</p>

# Quick start guide
Please check the package requirements in the [Installation Guide](https://ersilia.gitbook.io/ersilia-book/quick-start/installation). The next steps are a quickstart guide to installing Ersilia.

1. Create a conda environment and activate it
```
conda create -n ersilia python=3.7
conda activate ersilia
```
2. Clone this repository and install with pip
```
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
pip install -e .
``` 
3. Once the Ersilia Model Hub is installed, you can use the CLI to run predictions. First, select a model from the [Ersilia Model Hub](https://ersilia.io/model-hub/) and **fetch** it:
```
ersilia fetch chemprop-antibiotic
```
4. Generate a few (5) example molecules, to be used as input. The **example** command will generate the adequate input for the model in use
```
ersilia example chemprop-antibiotic -n 5 -f my_molecules.csv
```
5. Then, **serve** your model:
```
ersilia serve chemprop-antibiotic
```
6. And run the prediction **API**:
```
ersilia api -i my_molecules.csv -o my_predictions.csv
```
7. Finally, **close** the service when you are done.
```
ersilia close
```

Please see the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailed explanations.
<p align="right">(<a href="#top">back to top</a>)</p>

# Contribute
The Ersilia Model Hub is a Free, Open Source Software and we highly value new contributors. There are several ways in which you can contribute to the project:
* A good place to start is checking open [issues](https://github.com/ersilia-os/ersilia/issues). 
* If you have identified a bug in the code, please open a new issue using the bug template.
* If you want to incorporate a new model in the platform, open a new issue using the model request template or contact us using the following [form](https://www.ersilia.io/request-model)
* Share any feedback with the community using [GitHub Discussions](https://github.com/ersilia-os/ersilia/discussions) for the project
* Check our [Contributing Guide](https://github.com/ersilia-os/ersilia/blob/master/CONTRIBUTING.md) for more details

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) code of conduct.
<p align="right">(<a href="#top">back to top</a>)</p>

# License and citation
This repository is open-sourced under the [GPL-3 License](https://github.com/ersilia-os/ersilia/blob/master/LICENSE).
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.
<p align="right">(<a href="#top">back to top</a>)</p>

# About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.

[Help us](https://www.ersilia.io/donate) achieve our mission or [volunteer](https://www.ersilia.io/volunteer) with us!

<p align="right">(<a href="#top">back to top</a>)</p>
