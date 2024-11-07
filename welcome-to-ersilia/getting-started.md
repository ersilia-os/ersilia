---
description: >-
  Learn about the work of Ersilia and where to start using/contributing to our
  tools!
---

# Getting started

Ersilia develops and implements AI/ML tools for infectious disease research. This documentation will be useful if you are...

* A scientist looking to use some of our AI/ML platforms for your projects
* An open-source developer aiming to contribute to our tools
* A scientist developing AI/ML tools and wishing to incorporate them in our platform
* An Ersilia enthusiast looking forward to learning more about our work!

All of our work is openly available through our GitHub organisation page. Below we will summarise the main repositories and where to find the most important software tools. For a complete catalog of all our organisation repositories, please see [this table](https://airtable.com/app1iYv78K6xbHkmL/shrdboxOgKjy3vipl).

## Ersilia Model Hub

The Ersilia Model Hub is our main platform. It serves ready-to-use AI models related to the drug discovery cascade. Models can be browsed in our [website](https://ersilia.io/model-hub), and can run locally (see [Installation ](../ersilia-model-hub/installation.md)instructions) and we also offer a selection of them for online inference (please select those available Online through our [website](https://ersilia.io/model-hub)) as well as an [open service](https://github.com/ersilia-self-service) based on GitHub.

Detailed information about the Hub, its components and how to use it and contribute to its backend as well as contribute models can be found in the Ersilia Model Hub GitBook [section](broken-reference).

The repositories linked to the Hub are:

* [ersilia](https://github.com/ersilia-os/ersilia): platform backend.
* [ersilia-self-service](https://github.com/ersilia-os/ersilia-self-service): GitHub Action-based online inference for all models. Data and results are available publicly through GitHub issues. Please do not submit IP-sensitive data.
* [ersilia-assistant](https://github.com/ersilia-os/ersilia-assistant): LLM-based interface to easily interact with the Ersilia Model Hub.
* [ersilia-stats](https://github.com/ersilia-os/ersilia-stats): collection of statistics around the Hub and its usage, such as scientific publications, disease areas covered, etc.
* [eos-template](https://github.com/ersilia-os/eos-template): template for new model incorporation.
* [ersilia-pack](https://github.com/ersilia-os/ersilia-pack): model packaging for serving through FastAPI
* [ersilia-maintenance](https://github.com/ersilia-os/ersilia-maintenance): GitHub Action-based repository to check for integrity of the models within Ersilia
* [ersilia-maintained-inputs](https://github.com/ersilia-os/ersilia-maintained-inputs): standardised inputs for model testing
* [isaura](https://github.com/ersilia-os/isaura): data-lake to store predictions locally and in the cloud, saving computational time and resources
* [model-inference-pipeline](https://github.com/ersilia-os/model-inference-pipeline): isaura-based pipeline to store model inference results online, creating an open database of pre-calculated bioactivity and ADME properties of small molecules.
* eos repositories: repositories labelled as eosxxxx contain individual models. A full list of models, its identifier and relevant information is available in [this table](https://airtable.com/appgxpCzCDNyGjWc8/shrNc3sTtTA3QeEZu).

## ZairaChem

ZairaChem is an automated pipeline for ML model training. Read more about it in its dedicated [GitBook section](../chemistry-tools/automated-activity-prediction-models.md), as well as the [associated publication](https://www.nature.com/articles/s41467-023-41512-2) and [code repository](https://github.com/ersilia-os/zaira-chem). Coupled to ZairaChem, we have developed [Olinda](https://github.com/olinda), a model distillation framework to convert the high-performant, heavy ZairaChem models into portable .onnx models amenable for large library screening and online deployment.

## ChemSampler

ChemSampler is a pipeline based on the generative AI models available in the Ersilia Model Hub. Given a starting molecule, it performs several rounds of generative chemistry and produces a list of molecular candidates. ChemSampler can be constrained using several parameters. Please read its dedicated [GitBook section](../chemistry-tools/sampling-the-chemical-space.md) or check the code repository

## Other tools

As we work in a domain where IP-sensitive data is a common issue, we have developed a framework to ensure user queries are encrypted and thus that our models can be used by anyone. Chemxor is developed but not yet implemented in our ecosystem. Find more information about it in its dedicated [GitBook section](../privacy-preserving-ai/page-2.md) and associated [code repository.](https://github.com/ersilia-os/chemxor)

## Workshops and courses

As part of our mission we provide training in AI and Data Science to researchers across the Global South. All our trainings are documented and freely available. Check out the [Training Materials](broken-reference) section if you are interested, and have a look at the following code repositories:

* AI2050 courses: 2h introduction to Drug Discovery (code) and full week course for more advanced students (code), developed in collaboration with the H3D Foundation.
* Event Fund

## Paper-associated analyses and apps

*
* sars-cov-2-chemical-space
* osm...
* ligand-discovery...
