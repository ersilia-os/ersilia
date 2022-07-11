---
description: This tutorial explains how to incorporate models in the Ersilia Model Hub
---

# Incorporate models

{% hint style="danger" %}
This page is **work in progress**!
{% endhint %}

## Anatomy of the Ersilia Model Template

Each model in the Ersilia Model Hub is contained within an individual GitHub repository. The **Ersilia Model Template** in the [eos-template](https://github.com/ersilia-os/eos-template) repository is stored as a GitHub Template, so you can [create a new repository based on this template](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template).

Below, we describe the template files in detail. Note that we only explain the files that you need to modify; other files, like `pack.py`, do not need any modification from the model contributor.

### The `eos` identifier

Each model in the Ersilia Model Hub has an Ersilia Open Source (EOS) identifier. This identifier determines the name of the GitHub repository containing the model:

```
https://github.com/ersilia-os/[EOS_IDENTIFIER]
```

The `eos` identifier follows this regular expression: `"eos[1-9][a-z0-9]{3}"`. That is:

* The `eos` prefix, plus...
* one digit  (`1-9`) (the `0` is reserved for test models), plus...
* three alphanumeric (`a-z` and `0-9`) characters.

{% hint style="success" %}
A list of `eos` identifiers needs is provided in the [Ersilia Model Hub Spreadsheet](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=drive\_web\&ouid=114775674178390159004). You can read about this Spreadsheet [here](model-selection.md).
{% endhint %}

### The [`README`](https://github.com/ersilia-os/eos-template/blob/main/README.md) file

The `README.md` file is where the basic documentation about the model is provided. It must include the following fields:

**Title:** a self-descriptive model title (less than 70 characters)

**Model Identifiers:** a set of codes that identify the model in the Ersilia CLI:

* Ersilia Identifier (EOS ID): an automatically generated code (`eos` +  4 alphanumeric characters). Use the assigned identifier in the [Ersilia Model Hub Spreadsheet](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=sharing).
* Slug: a one-word or two-word (linked by a hypen) reference that acts as substitute to the EOSID
* Tags: labels to facilitate model search. For example, a model that predicts activity against malaria could have _P.Falciparum_ as tag. Select three relevant ones.

**Description**: minimum information about model type, results and the training dataset.

* Input: data format required by the model. Most chemistry related models, for example, will require molecules in SMILES format. If other input types, such as InChIKeys or peptide sequences are accepted, specify them.
* Output: unit and description of the model result. It is essential to specify if the model gives back a probability or a specific measure. F_or example, IC50_
* Model type: regression or classification
* Training set: number of compounds and link to the training dataset if available
* Mode of training: pretrained (the checkpoints where downloaded directly from a third party) retrained (the model was trained again using the same or a new dataset (please specify)) new model (if the model has been developed from scratch by Ersilia's contributors)

Results interpretation: provide a brief description of how to interpret the model results. _For example, in the case of a binary classification model for antimalarial activity based on experimental IC50, indicate the experimental settings (time of incubation, strain of parasite...) and the selected cut-off for the classification._

**Source Code:** this section must contain **all** relevant information about the model original authors, including a link to the publication if the model has been published in a peer reviewed journal or is in a preprint repositories, a link to the source code (typically, GitHub or GitLab) and a link to the model checkpoints directly, when available.

**License:** in addition to the `LICENSE.md` file, it is good practice to specify the Licenses in the README.md. All models in the Ersilia Model Hub are licensed under an open source license. Please make sure to abide by requirements of the original license when re-licensing or sub-licensing third-party author code. You can read more about how we deal with OS Licenses [here](https://ersilia.gitbook.io/ersilia-book/contributors/open-source-licences).

**History:** a short, numbered explanation of the modifications made in the source code, including the date of download and any steps taken to embed the model within the Ersilia Model Hub infrastructure.&#x20;

**About us:** all Ersilia repositories contain a final _About_ section, please keep the predefined version.

### The [`LICENSE`](https://github.com/ersilia-os/eos-template/blob/main/LICENSE) file

All original code written in contribution to Ersilia should be licensed under a GPLv3 License. The main License file of the repository, therefore, will be a GPLv3 as specified by [GitHub.](https://docs.github.com/en/communities/setting-up-your-project-for-healthy-contributions/adding-a-license-to-a-repository)

In addition, the license notices for code developed by third parties must be kept in the respective folders where the third-party code is found. Include an explanation in the `README` file, for example:  _The GPL-v3 license applies to all parts of the repository that are not externally maintained libraries. This repository uses the externally maintained library "Chemprop", located at /model and licensed under an_ [_MIT License_](https://github.com/ersilia-os/eos4e40/blob/main/model/LICENSE.md)_._

### The [`Dockerfile`](https://github.com/ersilia-os/eos-template/blob/main/Dockerfile) file

This file contains the installation instructions to install the model. Ersilia uses a `Dockerfile` to specify installation instructions since Docker contains the maximum level of isolation for the model. However, in most practical scenarios, a Docker image will not be necessary and a Conda environment, or even a Virtualenv environment, will suffice. The Ersilia CLI will decide which isolation level to provide depending on the content of the `Dockerfile`:

* If only `pip install` commands are specified, Virtualenv will be used.
* If only `pip install` and `conda install` commands are specified, Conda will be used.
* If other commands are specified (e.g. `sudo apt-get`), Docker will be used.

The `Dockerfile` available in the Ersilia Model Template is as follows:

```docker
FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit=2020.03

WORKDIR /repo
COPY ./repo
```

In this case, a Conda environment will be preferentially used to isolate the model. The first line of the `Dockerfile` indicates that this Conda environment will have **BentoML 0.11.0** installed on **Python 3.7**.

The Dockerfile can contain as many `RUN` commands as necessary, between the `MAINTAINER` and the `WORKDIR` lines.

### The [`model`](https://github.com/ersilia-os/eos-template/tree/main/model) folder

XX

### The [`.gitattributes`](https://github.com/ersilia-os/eos-template/blob/main/.gitattributes) file

We use Git LFS to store large files (over 100 MB), typically corresponding to model parameters. Files to be stored in Git LFS should be specified in the `.gitattributes` file. The current file will store in Git LFS all files in `csv`, `h5`, `joblib`, `pkl`, `pt` and `tsv` format.

```
*.csv filter=lfs diff=lfs merge=lfs -text
*.h5 filter=lfs diff=lfs merge=lfs -text
*.joblib filter=lfs diff=lfs merge=lfs -text
*.pkl filter=lfs diff=lfs merge=lfs -text
*.pt filter=lfs diff=lfs merge=lfs -text
*.tsv filter=lfs diff=lfs merge=lfs -text
```

You have to edit the `.gitattributes` file to ensure that all large files in your model are stored in Git LFS.

### The [`service`](https://github.com/ersilia-os/eos-template/blob/main/src/service.py) file

The service file is located in `src/service.py`. It contains the necessary code to facilitate model bundling with BentoML.

There are three main classes in the `service` file, namely `Model`, `Artifact` and `Service`.

#### The `Model` class

This class is simply a wrapper for the AI/ML model. The most important method of the `Model` class is the `predict` method.

{% hint style="info" %}
You will have to rename the `predict` method into something else if your model does not predict. For example, for many models the method can be renamed as `calculate.`&#x20;
{% endhint %}

{% hint style="success" %}
Multiple methods are allowed. For example, a model may have a `predict` and an e`xplain` method.
{% endhint %}

#### The `Artifact` class

XX

#### The `Service` class

XX

## Steps for model incorporation

In this tutorial, we will follow the example of a very simple model, related synthetic accessibility scoring.

#### Find model code and parameters

XX

#### Run the code _outside_ Ersilia

XX







