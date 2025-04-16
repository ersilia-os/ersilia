---
description: >-
  Learn about the work of Ersilia and where to start using/contributing to our
  tools
---

# Ersilia's ecosystem

Ersilia develops and implements AI/ML tools for infectious disease research. This documentation will be useful if you are...

* A **chemist or biologist** looking to use some of our AI/ML platforms for your projects.
* An open-source **developer** aiming to contribute to our tools.
* A **data scientist** developing AI/ML tools and wishing to incorporate them in our platform.
* An Ersilia **enthusiast** looking forward to learning more about our work.

All of our work is openly available through our [GitHub organisation page](https://github.com/ersilia-os). Below we will summarise the main repositories and where to find the most important software tools. For a complete catalog of all our organisation repositories, please see [this table](https://airtable.com/app1iYv78K6xbHkmL/shrdboxOgKjy3vipl).

## Tool repositories

### The Ersilia Model Hub

The [Ersilia Model Hub](https://ersilia.io) is our main platform. It serves ready-to-use AI models related to the drug discovery cascade. Models can be browsed in our [website](https://ersilia.io/model-hub), and can run locally (see [Installation ](broken-reference)instructions) and we also offer a selection of them for online inference (please select those available Online through our [website](https://ersilia.io/model-hub)) as well as an [open service](https://github.com/ersilia-self-service) based on GitHub.

Detailed information about the Ersilia Model Hub, its components and how to use it and contribute to its backend as well as contribute models can be found in [this section](broken-reference). Developers may look into the [API documentation](https://ersilia-os.github.io/ersilia) for an in-depth view of the code.

The repositories linked to the Ersilia Model Hub are:

* [ersilia](https://github.com/ersilia-os/ersilia): This is the main repository, corresponding to a CLI to fetch and run models locally.
* [eos-template](https://github.com/ersilia-os/eos-template): Template for new model incorporation. This template uses GitHub Actions workflows specified in [ersilia-model-workflows](https://github.com/ersilia-os/ersilia-model-workflows).
* [ersilia-stats](https://github.com/ersilia-os/ersilia-stats): Collection of statistics around the Hub and its usage, such as scientific publications, disease areas covered, etc.
* [ersilia-maintenance](https://github.com/ersilia-os/ersilia-maintenance): GitHub Actions-based repository to check for integrity of the models within Ersilia.
* [ersilia-self-service](https://github.com/ersilia-os/ersilia-self-service): GitHub Action-based online inference for all models. Data and results are available publicly through GitHub issues. Please do not submit IP-sensitive data.
* [ersilia-assistant](https://github.com/ersilia-os/ersilia-assistant): LLM-based interface to easily interact with the Ersilia Model Hub.
* [ersilia-pack](https://github.com/ersilia-os/ersilia-pack): Model packaging for serving through FastAPI.
* [ersilia-maintained-inputs](https://github.com/ersilia-os/ersilia-maintained-inputs): Standardised inputs for model testing.
* [model-inference-pipeline](https://github.com/ersilia-os/model-inference-pipeline): Pipeline to store model inference results in AWS, creating an open database of pre-calculations (cache).
* [eos repositories](https://github.com/ersilia-os/eos4e40): repositories labelled with an Ersilia (eos) identifier contain individual models. A full list of models, their identifiers and relevant information is available in [this table](https://airtable.com/appgxpCzCDNyGjWc8/shrNc3sTtTA3QeEZu).

### ZairaChem

ZairaChem is an automated pipeline for ML model training. Read more about it in its dedicated [section](../chemistry-tools/automated-activity-prediction-models.md), as well as the [associated publication](https://www.nature.com/articles/s41467-023-41512-2) and code repository ([zaira-chem](https://github.com/ersilia-os/zaira-chem)). Coupled to ZairaChem, we have developed Olinda, a model distillation framework to convert the high-performant, heavy ZairaChem models into portable ONNX models amenable for large-scale calculations and online deployment ([olinda](https://github.com/ersilia-os/olinda)).

### ChemSampler

ChemSampler is a pipeline based on the generative AI models available in the Ersilia Model Hub. Given a starting molecule, it performs several rounds of generative chemistry and produces a list of molecular candidates. ChemSampler can be constrained using several parameters. Please read its dedicated [GitBook section](../chemistry-tools/sampling-the-chemical-space.md) or check the code repository ([chem-sampler](https://github.com/ersilia-os/chem-sampler)).

{% hint style="danger" %}
ChemSampler is still under development. Please open a [GitHub Issue](https://github.com/ersilia-os/chem-sampler/issues) if you want to use this tool, and we will try to assist you accordingly.&#x20;
{% endhint %}

## Workshops and courses

As part of our mission we provide training in AI and Data Science to researchers across the Global South. All our trainings are documented and freely available. Check out the [Training Materials](broken-reference) section if you are interested, and have a look at the following code repositories:

* AI2050 courses: 2h introduction to Drug Discovery ([ai2050-h3d-symposium-workshop](https://github.com/ersilia-os/ai2050-h3d-symposium-workshop)) and full week course for more advanced students ([ai2050-dd-workshop](https://github.com/ersilia-os/ai2050-dd-workshop)), developed in collaboration with the H3D Foundation.
* Event Fund: A one-week course we developed in collaboration with the H3d Centre and the support of the Wellcome Trust and Code for Science and Society.
* Python 101: An introduction to Python programming language geared to scientists (focusing on data analysis, plotting and basic pythonic operations; [python101](https://github.com/ersilia-os/python101)). Inspired by the Carpentries!

## Research-associated analyses

In addition to our software tools, we have a number of repositories related to scientific research projects. Those repositories typically contain the necessary data and code to reproduce an analysis reported in a research paper. For a full overview of our research projects and publications please have a look at our [website](https://ersilia.io/publications). Below are a few exemplary projects, finalized or in current development:

* ADDA4TB: Targeted protein degradation for Mycobacterium tuberculosis, in collaboration with Stellenbosch University ([mtb-targeted-protein-degradation](https://github.com/ersilia-os/mtb-targeted-protein-degradation)).
* GRADIENT Pharmacogenetics in Africa: Analysis of potential pharmacogenes related to antimalaria and anti-TB drugs, in collaboration with H3D ([pharmacogx-embeddings](../chemistry-tools/automated-activity-prediction-models/accurate-automl-with-zairachem.md)).
* SARS-CoV-2 Chemical Space: Analysis of the chemical space associated with curated COVID-19 therapy data, done in collaboration with UB-CeDD ([sars-cov-2-chemical-space](https://github.com/ersilia-os/sars-cov-2-chemspace)).
