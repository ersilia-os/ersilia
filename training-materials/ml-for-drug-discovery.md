---
description: >-
  This page describes the training materials available at Ersilia for scientists
  starting in the field of Artificial Intelligence and Machine Learning for Drug
  Discovery
---

# ML for Drug Discovery

## Course Overview

The goal of this course is to provide a first insight into the use of AI/ML tools for drug discovery, with a strong focus on the prediction of bioactivities against pathogens. At course completion, students will be able to:

* Read and understand publications in the domain of computational biology
* Know where to look for public datasets to build AI/ML models
* Leverage open source tools like GitHub
* Expand on basic concepts of Python programming and Jupyter Notebooks
* Distinguish between several cloud and local-based computing systems
* Use already existing AI/ML models and apply them to their research

The course is geared towards graduate students and postdoctoral researchers with a background in biology, biomedicine or chemistry who want to focus their work in the exciting field of computational biology and data science! Preferably course participants should have active ongoing research projects.

{% hint style="warning" %}
This course is not aimed at computer scientists savvy in programming and data management. We provide foundations to understand how to apply AI/ML to research projects, not a deep-dive into AI/ML development.
{% endhint %}

## Contents

The course is organised in 5 modules. Each module focuses on one aspect of drug discovery and AI/ML, and showcases a specific set of tools. Surveys are used throughout the workshop to encourage participation and help recap the important concepts.

The course duration is 10 full days.

### Module 0: Introduction to drug discovery

Basic steps of the drug discovery cascade and how can AI/ML models aid in each one of them. Module 0 also introduces at a general level the steps needed to build an AI/ML model, which will be expanded in the following modules.

* Duration: 1 day
* Material:&#x20;
  * [Presentation](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/ersilia-intro-workshop-m0.pdf)
  * Knowledge surveys [0a](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m0a.pdf) and [0b](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m0b.pdf)

### Module 1: Using AI models for drug discovery

This model brings a hands-on exercise where students are faced with a real-life problem. We provide two lists of compounds, one extracted from ChEMBL and one extracted from Coconut, and the participants must select the best compounds for experimental antimalarial screening. To do so, they must use the models in the Ersilia Model Hub, select a few relevant ones and run predictions using the Graphical User Interface of Ersilia (advanced CLI use will be presented in M2 and M4)

* Duration: 2 days
* Material:
  * [Presentation](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/ersilia-intro-workshop-m1.pdf)
  * Knowledge surveys [1a](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m1a.pdf) and [1b](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m1b.pdf)
  * Antimalarial screening [notebook](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/notebooks/m1\_antimalarial\_exercise.ipynb)

### Module 2: Setting up the computational environment

The student learns about tools for data science, including basic Python and bash programming, Git and GitHub as well as cloud environments like Google Colab and GitHub Codespaces. The choice of cloud environment will depend on the capabilities (computer equipment, internet speed...) at the site of course delivery. The extension allows for deepening the study of Python language.

* Duration: 1 day + 1 day extension
* Material:
  * [Presentation](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/ersilia-intro-workshop-m2.pdf)
  * Knowledge surveys [2a](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m2a.pdf), [2b](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m2b.pdf) and [2c](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m2c.pdf)
  * Python extension [notebooks](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/notebooks/m2\_python.ipynb) and [exercise notebooks](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/notebooks/m2\_python\_exercises.ipynb)

### Module 3: Training an AI model for bioactivity prediction

This module starts with the presentations by the students of their own research projects. The presentations must be simple and try to identify a few aspects of the project where AI/ML modelling could help, based on the learnings from the previous modules. Following this introduction, the course facilitators will identify a few datasets relevant to the area of study of the students and prepare them for building an AI model. An example is given with a template notebook ready to run on Colab. If time allows, course facilitators can expand the work to regression models, otherwise the content will be focused on classifications.

* Duration: 3 days
* Material:
  * [Presentation](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/ersilia-intro-workshop-m3.pdf)
  * [Exercise example](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/ersilia-intro-workshop-m3-exercises.pdf)
  * Knowledge survey [3](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/surveys/ersilia-intro-m3.pdf)
  * Template [notebook](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/notebooks/m3\_building\_a\_model.ipynb) for AI model training

### Module 4: The Ersilia Model Hub

The course concludes with the addition of the models built in the class to the Ersilia Model Hub. The students learn how to install and run Ersilia in their own systems, ideally in GitHub Codespaces first, as the GitHub repository is set up already for running Ersilia.

* Duration: 2 days
* Material:
  * [Presentation](https://github.com/ersilia-os/ersilia-intro-workshop/blob/main/presentations/ersilia-intro-workshop-m4.pdf)
  * [Ersilia documentation](../ersilia-model-hub/antibiotic-activity-prediction.md)

{% hint style="info" %}
If time allows, the best to finish the course is a presentation in small groups of the work developed during M3 and M4
{% endhint %}

## Resources

All the course material and content can be found in the [Ersilia Introductory Workshop](https://github.com/ersilia-os/ersilia-intro-workshop) GitHub repository.

If new models are developed as part of the module 3 practice, those will be also added to the GitHub repository. The course surveys are developed using Mentimeter, but another software can be used for facilitation if preferred.

#### Additional courses

We have developed a 4-day training course on the applications of AI/ML to infectious disease drug discovery aimed at chemists and biologists who want to apply the insights of AI/ML models to their work.&#x20;

The course contents are divided as follows:

* Introduction to Drug Discovery and working with Chemistry data
* Introduction to ML basics and model outputs
* ADMET tools
* Open Science and generative models

Each block content is organised as a full day session, including keynote lectures, hands-on training and breakout sessions. If you cannot attend a live course, we recommend using the materials for self-learning. All relevant content and links, as well as extra materials, can be found in the [Event Fund](https://ersilia.gitbook.io/event-fund) webpage.

All content is released under a CC-BY-4 license.
