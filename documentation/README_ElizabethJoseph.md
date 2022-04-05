
# Ersilia Model Hub


[![donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


[![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Codestyle](https://img.shields.io/badge/code%20style-black-000000.svg?logo=Python&logoColor=white)](https://github.com/psf/black) [![Gitpod](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia)


![Ersilia Logo](https://softr-prod.imgix.net/applications/5ce288d5-4600-42df-9d5f-8617b023a3e3/assets/f7b7b7e7-8f53-425e-8504-f722cfa0f503.png)

## Table of content
1. [Project Description](https://github.com/ersilia-os/ersilia#project-description)
2. [Getting started](https://github.com/ersilia-os/ersilia#getting-started)
3. [Contribution](https://github.com/ersilia-os/ersilia#contribution)
4. [Language](https://github.com/ersilia-os/ersilia#language)
4. [Roadmap](https://github.com/ersilia-os/ersilia#roadmap)
5. [License and citation](https://github.com/ersilia-os/ersilia#license-and-citation)
6. [About us](https://github.com/ersilia-os/ersilia#about-us)
8. [Follow us](https://github.com/ersilia-os/ersilia#Follow-us)  


The Ersilia Model Hub is the main project of the [Ersilia Open Source Initiative](https://ersilia.io/). The aim is to provide a platform for a user-friendly deployment of AI/ML models, where scientists can browse through the assets, identify those which is relevant to their research and obtain predictions without the need of coding expertise. Currently, most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who does not have the necessary expertise. This gap becomes even larger in Low and Middle Income Country institutions where access to bioinformatic facilities or data science experts are scarce. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.  


The platform is aimed at helping researchers identify drug candidates for orphan and neglected diseases, design molecules from the beginning and understand the mechanisms of action or anticipate adverse side effects. The ultimate goal of Ersilia is to lower the barrier to drug discovery, encouraging academic groups and companies to pursue the development of new medicines following the principles of Open Science.  

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgement) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.
### A schematic view of the Ersilia Model Hub.

![Ersilia Model Hub Image](https://2591732297-files.gitbook.io/~/files/v0/b/gitbook-legacy-files/o/assets%2F-Mj44wxA7bU1hQH19m8I%2F-MjPA_JIOzuKGPtPDw1U%2F-MjPDmAAJkYPLJmKKSD1%2FErsilia_Hub-02.png?alt=media&token=8a876edc-c02e-400c-80c5-b6f3c4060c21)

* Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)
# Getting started
## Installation
There are a few third-party **pre-requisites** and you need to have them installed on your computer before proceeding with the installation of the Ersilia tool. Below is a list of the required programs, you'll be redirected to the websites for download.
* [Python 3.7 and above](https://www.python.org/)
* [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) / [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
* [Docker](https://www.docker.com/)
* [GitHub CLI](https://cli.github.com/manual/installation)
* [Git LFS](https://git-lfs.github.com/)



## Installing on Linus and Macosx
1. Open the Conda terminal or environment  
 **Create a Conda environment**  
    
```
conda create -n ersilia python=3.7
```

**Activate the environment**
```
conda activate ersilia
```
2. Install the Ersilia Python package
**Clone from GitHub**  
```
git clone https://github.com/ersilia-os/ersilia.git
```
```
cd ersilia
```
**Install with pip (use -e for developer mode)**
```
pip install -e .
```
3. Quickly check that the CLI works on your terminal
**See ersilia CLI options**
```
ersilia --help
```

## The Isaura data lake
We highly recommend installation of the [Isaura](https://github.com/ersilia-os/isaura) data lake. With Isaura, you will be able to cache your model predictions (i.e. store them in your local computer). Isaura is a relatively light Python package:  

**Clone from github**  
```
git clone https://github.com/ersilia-os/isaura.git
cd isaura
```
**Install with pip (use -e for developer mode)**  
```
pip install -e .
```

## Installing on Windows 
We are not testing Windows installation consistently. If you encounter problems, please reach out to us at hello@ersilia.io.  
We recommend that you install Ubuntu on your Windows machine. This can be now done very easily with WSL. You will need at least Windows 10.  

**Open a Power Shell with Admin permissions and type:**
```
wsl --install
```
Then simply [install the Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) on Windows.  
Inside the Ubuntu terminal, you can now follow the installation instructions [above](https://ersilia.gitbook.io/ersilia-book/quick-start/installation#install-on-linux-and-macosx).

For better installation guide of the Ersilia API Model, follow the detailed steps on the [ Ersilia Model Book](https://ersilia.gitbook.io/ersilia-book/quick-start/installation/).

## Testing Ersilia Model Hub
The Ersilia Model Hub can be used to test a list of computational assets that you can use, for instance, to calculate or predict properties for your molecules of interest.
To test our model, we will be calculating **Molecular Weight.**

First, you have to **fetch** the model from our repository. The molecular weight calculator is named ```molecular-weight:```  

**Download and Install**
```  
ersilia fetch molecular-weight
```
Then, you can serve this tool so that it becomes available as an API:  
**Start a service**  
```
ersilia serve molecular-weight
```
Our **API** of interest is ```calculate.``` We can use it to calculate the molecular weight of Aspirin.    

**Calculate molecular weight of Aspirin using SMILES**
```
ersilia api calculate -i "CC(=O)OC1=CC=CC=C1C(=O)O"
```
Now that we are done, we can **close** the service.  
  
**Close service**
```
ersilia close
```
Please see the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/) for more examples and detailed explanations.
# Contributions
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! Below are the ways in which you can contribute to the development of this project:
* To improve our tool as a developer, check out the [Issues](https://github.com/ersilia-os/ersilia/issues)
* To include a Model you developed in the Hub, kindly [Contact us](https://ersilia.io/) or open an [Issue](https://github.com/ersilia-os/ersilia/issues)
* If you are a scientist with a cool dataset, also [Contact us](https://ersilia.io/) or open an [Issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a Model based on the data!
* If there is a third-party model you have identified and would like to see it in the Hub, open an [Issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.
* If you have any questions or concerns about this Model, kindly [Contact us](https://ersilia.io/)
The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.


Many thanks to these contributors for their work.
[![](https://github.com/GemmaTuron.png?size=50)](https://github.com/GemmaTuron) **Gemma Turon**
[![](https://github.com/miquelduranfrigola.png?size=50)](https://github.com/miquelduranfrigola) **Miquel Duranfrigola**


*Contributions of any kind are welcomed.*
# Language
The Languages used in the Ersila Model Hub are:
* **Python**
* **Dockerfile**
# Roadmap
We are working to grow the Hub organically and responding to our users needs. Here a detail of the next features to come, stay tuned!

1. Deployment for Windows System (expected: February 2022)
2. Automated third-party model contributions (expected: March 2022)
3. Possibility to run lite models online (expected: May 2022)
# License and Certification
This repository is open-sourced under the MIT License. Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.
# About us
The [Ersilia Open Source Initiative](https://ersilia.io/) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

You can support us via [Open Collective](https://github.com/opencollective.com/ersilia).
# Follow us
* Join the Ersilia [Slack Workspace](https://ersilia-workspace.slack.com/)  
* Subscribe to our [mailing list](https://twitter.com/ersiliaio)
