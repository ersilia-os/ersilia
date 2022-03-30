# Welcome to the Ersilia Model Hub!

[![Donate](https://img.shields.io/badge/Donate-PayPal-green.svg)](https://www.paypal.com/uk/fundraiser/charity/4145012) [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](code_of_conduct.md) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![documentation](https://img.shields.io/badge/-Documentation-purple?logo=read-the-docs&logoColor=white)](https://ersilia.gitbook.io/ersilia-book/) [![EOSI](https://circleci.com/gh/ersilia-os/ersilia.svg?style=svg)](https://circleci.com/gh/ersilia-os/ersilia) [![PyPI version fury.io](https://badge.fury.io/py/ersilia.svg)](https://pypi.python.org/pypi/ersilia/) [![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/) [![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/ersilia-os/ersilia) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity) [![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg) [![macOS](https://svgshare.com/i/ZjP.svg)](https://svgshare.com/i/ZjP.svg) [![Windows](https://svgshare.com/i/ZhY.svg)](https://svgshare.com/i/ZhY.svg)

![logo](https://raw.githubusercontent.com/ersilia-os/ersilia/master/assets/Ersilia_Plum.png)

### Table of Contents:
1. [Project Description](#project-description)
2. [Installation Guide](#installation-guide)
3. [Contribute](#contribute)
4. [Roadmap](#roadmap)
5. [License and citation](#license-and-citation)
6. [About us](#about-us)
7. [FAQs](#faqs)
8. [Related Projects](#related-projects)
## Project Description
The Ersilia Model Hub is a patient zero of [Ersilia Open Source Initiative](https://ersilia.io). The aim is to provide a platform for a user-friendly deployment of AI/ML models, where scientists can browse through the assets, identify those which is relevant to their research and obtain predictions without the need of coding expertise. Currently, most of the tools developed, even when published and fully open-sourced, remain unusable by a large majority of the scientific community who does not have the necessary expertise. This gap becomes even larger in Low and Middle Income Country institutions where access to bioinformatic facilities or data science experts are scarce. With this project, we hope to democratize access to this expertise and support research into neglected and infectious diseases.

The models embedded in the hub include both models published in the literature (with appropriate third party acknowledgment) and models developed by the Ersilia team or contributors. All assets are open source, but please do check the Licenses of each model before using them.

* Read more about the project better at the [Ersilia Book](https://ersilia.gitbook.io/ersilia-book/)
* Available models can be checked at [Ersilia Model Hub](https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP)

## Installation Guide
### Prerequisites
- [Docker](https://docs.docker.com/engine/install/)
- [Python](https://www.python.org/downloads/)
- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
- [Git LFS](https://git-lfs.github.com/)
- [Git CLI](https://cli.github.com/manual/installation)
### Install on Linux and MacOSX
**Step 1: Set up conda environment**
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
#### Set up [Isaura](https://github.com/ersilia-os/isaura)
Clone from Github
```
git clone https://github.com/ersilia-os/isaura.git
cd isaura
```
### Usage
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
## contribute
A small team of Ersilia staff and volunteers develops and maintains the Ersilia Model Hub, and any help is greatly appreciated! There are a number of ways you may help with the project: 

1 Create a GitHub account (if you don't have one).

2 Go through documentation on basic knowledge on how to use, write and contribute on GitHub , especially how to make a pull request and how to commit changes. You can find the video link below very helpful.
https://youtu.be/PQsJR8ci3J0.

3 Go through the Ersilia's code of conduct and familiarize yourself with it. Please do ensure you adhere to the code of conduct during your contribution phase.

4  When contributing to a project, the first step is to go to the project site and pick an issue that you'd want to work on and that fits your skill set. Before tackling larger and more difficult issues as a first-time contributor, you should try to find something modest and reasonably simple to utilize as a pleasant entry into the project. On your first few contributions, don't try to go too deep!

5 If you see an issue you are interested in tackling, you can indicate interest in the comment section to be assigned that issue. You can also open issues you feel should be added.

6  After making a contribution ensure you record your contribution on the outreachy website. https://www.outreachy.org/apply/project-selection/

7  If you have developed a model and would like to include it in the Hub to increase its visibility and usability, please https://ersilia.io/ or open an issue. We are currently working on an automated contribution tool to facilitate the process.

8 If you are a scientist with a cool dataset, also contact us or open an issue as we might be interested in developing a model based on the data!

9  If there is a third-party model you have identified and would like to see it in the Hub, open an issue with the relevant information and we will get back to you as soon as possible.

The Ersilia Open Source Initiative adheres to the [Contributor Covenant](https://ersilia.gitbook.io/ersilia-wiki/code-of-conduct) guidelines.

## Roadmap
We are working to grow the Hub organically and responding to our users needs. Here a detail of the next features to come, stay tuned!

- [x] Deployment for Windows System (expected: February 2022)
- [x] Automated third-party model contributions (expected: March 2022)
- [ ] Possibility to run lite models online (expected: May 2022)
## License and citation
This repository is open-sourced under the MIT License.
Please [cite us](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) if you use it.
## About us
The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization incorporated with the Charity Commission for England and Wales (number 1192266). Our mission is to reduce the imbalance in biomedical research productivity between countries by supporting research in underfunded settings.

## FAQs
**Can I modify your code?**  
Definitely! Our code is released under an MIT license. You can copy, download and modify the code provided you maintain the license notice. Code from third-party authors is licensed according to the original license, in that case, please check carefully the license restrictions.

**Is there a subscription fee?**  
The Ersilia Model Hub is available for everyone to download and run in their computers. You can become a subscriptor and contribute to our mission via Open Collective.

**What is your legal structure?**  
The Ersilia Open Source Initiative is a Charitable Incorporated Organization recognized by the Charity Commission of England and Wales (number 1192266)

**Can I get a tax deduction on my donation?**  
Yes, if you are based in the UK your donations are Tax Deductible. Help us further by completing the GiftAid declaration and your gift will be matched at a rate of 25cts for each pound.

## Related Projects
- [Capacity Building](https://www.ersilia.io/capacity-building): Resources and courses tailor-made to the needs of our collaborators.
- [Automated modelling for low-resourced settings](https://github.com/ersilia-os/zaira-chem): End-to-end package to train and deploy AI/ML models based on chemistry data.
- [Ersilia Model Hub](https://ersilia.gitbook.io/ersilia-book/): An open-source repository of AI/ML models for global health and neglected diseases.
- [Privacy-preserving AI](https://www.ersilia.io/coming-soon): Development of an encryption tool to enable AI/ML model training using IP-sensitive data.
- [Virtual screening of African natural products](https://www.ersilia.io/coming-soon): Analysis of the natural product landscape to identify new anti infectives.
- [Open Source Malaria Series 4](https://github.com/ersilia-os/osm-series4-candidates-2): Design of new antimalarials using reinforcement learning algorithms.

