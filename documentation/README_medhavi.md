# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity) [![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg) [![macOS](https://svgshare.com/i/ZjP.svg)](https://svgshare.com/i/ZjP.svg) [![Windows](https://svgshare.com/i/ZhY.svg)](https://svgshare.com/i/ZhY.svg)

![logo](https://raw.githubusercontent.com/ersilia-os/ersilia/master/assets/Ersilia_Plum.png)

### Table of Contents:
1. [Project Description](#project-description)
2. [Installation](#installation)
3. [Contribution-guidelines](#contribution-guidelines)
4. [Roadmap](#roadmap)
5. [About us](#about-us)
6. [License](#license)

## Project Description

Ersilia is a nonprofit organisation that supports research for infectious and neglected diseases. Our goal is to facilitate adoption of artificial intelligence tools to accelerate the development of new medicines. We focus on establishing collaborations in low-resourced settings where the costs of drug discovery are prohibitive.
The Ersilia Model Hub provides pre-trained Artificial Intelligence models to scientists so that they can integrate these new technologies within their daily experimental tasks. FOSS software release should be accompanied by appropriate documentation, usage examples and community guidelines. In this project, we want to focus on working around the different documentation issues for the Ersilia Model Hub. This project is focused on developing the technical documentation for installation, usage and contribution and/or, depending on the intern’s interests, we will focus on strengthening the community guidelines, working on materials for outreach and dissemination (blogposts, scientific publications…).


This project is inspired by the participation of the Ersilia Open Source Initiative in three programs this year: The [Digital Infrastructure Incubator](https://github.com/matiassingers/awesome-readme), focused on FOSS community and governance, the [Open Life Sciences](https://openlifesci.org/), , for implementing Open Science practices at all levels and the [Software Sustainability Fellowship](https://www.software.ac.uk/programmes-and-events/fellowship-programme), oriented towards the sustainability of research software tools. The intern will gain access to all documentation and lessons learned during these programs and help implement them.

## Installation

### Install on Linux and MacOSX
**Step 1 : create a conda environment**
```
conda create -n ersilia python=3.7
```
**Step 2: Activate the environment**
```
conda activate ersilia
```
**Step 3: Install ersilia python package**
Clone from Github
```
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
```
Install with pip
```
# install with pip (use -e for developer mode)
pip install -e .
```
Test CLI    
```
ersilia --help
```
## Windows Installation
Installing Ersilia on a system running MacOS or a Linux is relatively straightforward because of UNIX-like shell and also how the application has been developed.
But that doesn't mean that running a Windows OS will exclude you from installing and using Ersilia.

For you to install Ersilia on your Windows machine, you need to enable the feature of Windows to install a support for Linux shell. This is called Windows Subsystems for Linux. 

Windows Subsystem for Linux (WSL) lets developers run a GNU/Linux environment -- including most command-line tools, utilities, and applications -- directly on Windows, unmodified, without the overhead of a traditional virtual machine or dual-boot setup. ([source](https://docs.microsoft.com/en-us/windows/wsl/))

For installation of WSL on older windows builds, you can look up the steps [here](https://docs.microsoft.com/en-us/windows/wsl/install-manual) or follow the steps below.

**1. Enable the Windows Subsystem for Linux**

Open PowerShell as Administrator (Start menu > PowerShell > right-click > Run as Administrator) and enter this command:
```
dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
```

**2. Enable Virtual Machine feature**
Open Powershell as Administrator just like the first step and enter the command below:
```
dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
```

**3. Restart your machine and continue from the next step**

**4. Download and Install the Linux kernel update package**
Download the package from [here](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi) if your machine is Windows x64.

Run the update package downloaded above. (Double-click to run - you will be prompted for elevated permissions, select ‘yes’ to approve this installation.)

 **5. Install the Linux distro of your choice**
You can install Ubuntu, or any other Linux distribution you prefer from the Microsoft Store

**6. Setup username and password for Linux installation**

**7. Run the following commands to get Linux installation ready to Ersilia installation**
```
sudo apt-get update
```
You'll provide the password your setup the Linux installation with whenever you run a command with sudo
```
sudo apt install build-essentials
```
**8. Install Miniconda**
```
mkdir -p ~/miniconda3
```
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
```
```
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
```

**9. Install Git**
```
sudo apt-get install git
```
**10. Create a virtual env with Miniconda**
```
conda create -n ersilia python=3.7
```
```
conda activate ersilia
```
**11. Clone and Install Ersilia**

Clone from github
```
git clone https://github.com/ersilia-os/ersilia.git
```

Change working directory to the repository from GitHub
```
cd ersilia
```

Install with pip (use -e for developer mode)
```
pip install -e .
```

**12. Final. Confirm your installation**
```
ersilia --help
```

**Set up [Isaura]**(https://github.com/ersilia-os/isaura)
Clone from Github
```
git clone https://github.com/ersilia-os/isaura.git
cd isaura
```
**Usage
**Fetch** a model
```
ersilia fetch <model name>
```
**Generate** example modules
```
ersilia example <model name> -n 5 -f my_molecules.csv
```
**Serve** your model
```
ersilia serve <model name>
```
**Run** Prediction API
```
ersilia api -i my_molecules.csv -o my_predictions.csv
```
**Close** the service when you are done
```
ersilia close

```
## Contribution-guidelines
The Ersilia Model Hub is developed and maintained by a small team of Ersilia employees and volunteers, and any contribution is highly valued! There are several ways in which you can contribute to the project:
- If you are a developer, check the [issues](https://github.com/ersilia-os/ersilia/issues) and help us to improve the tool
- If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues). We are currently working on an automated contribution tool to facilitate the process.
- If you are a scientist with a cool dataset, also [contact us](https://ersilia.io) or open an [issue](https://github.com/ersilia-os/ersilia/issues) as we might be interested in developing a model based on the data!
- If there is a third-party model you have identified and would like to see it in the Hub, open an [issue](https://github.com/ersilia-os/ersilia/issues) with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

## Roadmap
We are working to grow the Hub organically and responding to our users needs. Here a detail of the next features to come, stay tuned!

- [x] Deployment for Windows System (expected: February 2022)
- [x] Automated third-party model contributions (expected: March 2022)
- [ ] Possibility to run lite models online (expected: May 2022)

## About Us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

##  License
This repository is open-sourced under the MIT License.
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.
