---
description: >-
  This section shows how to download models and run predictions using the
  Ersilia Model Hub
---

# Model usage

You can explore the available models through our website or by running the following command in the CLI:

```
# display ready to use models
ersilia catalog
```

Each model is identified by:

* EOS-ID: `eos[1-9][a-z0-9]{3}`
* Slug: 1-3 word reference for the model

In this case example, we show how to run predictions based on the AI/ML model developed in the paper _Retrosynthetic accessibility score (RAscore) – rapid machine learned synthesizability classification from AI driven retrosynthetic planning_ by [Thakkar et al, 2021](http://dx.doi.org/10.1039/D0SC05401A). The RA score is particularly useful to pre-screen large libraries of compounds, for example those produced by generative models.

## Use model through CLI

### Fetch model and install it locally

The first step is to download the model to your local device and install it along with its dependencies. By default, a `~/eos` directory (for Ersilia Open Source) will be created in your `HOME`. This folder will contain all fetched models along with additional files to manage the AI/ML content available locally.

To download and install the RA Score prediction model, simply use the `fetch` command. In the Ersilia Model Hub, the RA Score prediction model has the **identifier** `eos2r5a`  and the **slug** `retrosynthetic-accessibility`. You can use either one to refer to this model all of the commands below

```bash
# fetch model from remote repository using slug ...
ersilia fetch retrosynthetic-accessibility
# ... or using ersilia identifier
ersilia fetch eos2r5a
```

### Get model information

Once the model is downloaded, you can get more information through the model card:

```bash
# display model card using slug...
ersilia card retrosynthetic-accessibility
# ... or using ersilia identifier
ersilia card eos2r5a
```

{% hint style="info" %}
We do our best to keep the user away from the [dependency hell](https://en.wikipedia.org/wiki/Dependency\_hell). Models are **automatically** installed with the necessary degree of isolation from the system. All models are available through GitHub and also as Docker Images
{% endhint %}

### Serve model

Once the model has been fetched, it should be ready to be used. A model in the Ersilia Model Hub can be thought of as a set of APIs. You can serve the model like this:

```bash
# serve model
ersilia serve retrosynthetic-accessibility
```

A URL will be prompted as well as a process id (PID). These can be relevant if you are an advanced user and want to have low-level control of the tool. The most important is, however, the list of **available APIs**. By default, all models use the `run` API.

### Make predictions

The RA Score prediction model takes **chemical structures** as input and provides a score (ranging from 0 to 1). The higher the score, the more synthetically accessible the molecule is predicted to be.

Ideally, in the chemistry models, the input molecules are specified as **SMILES** strings. SMILES strings can be easily found online. For instance, we can find an antibiotic, Halicin, in [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/Halicin#section=Canonical-SMILES), and then predict its retrosynthetic accessibility as follows:

```bash
# Halicin
ersilia api run -i "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]"
# also works with the run command directly
ersilia run -i "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]"
```

{% hint style="warning" %}
It is also possible to use [InChIKey](https://pubchem.ncbi.nlm.nih.gov/compound/Halicin#section=InChI-Key) or even molecule name (through the [Chemical Identifier Resolver](https://cactus.nci.nih.gov/chemical/structure)) instead of SMILES. Ersilia will take care of this automatically. However, please take into account that this requires an internet connection and will slow down the process, as requests to external tools are necessary.
{% endhint %}

You can make **multiple predictions** in batch mode. This is typically much faster than running predictions one by one in a loop. For instance, we can predict the RA Score of Halicin and [Ibuprofen](https://pubchem.ncbi.nlm.nih.gov/compound/Ibuprofen#section=Canonical-SMILES).

```bash
# Halicin and Ibuprofen
ersilia api run -i "['C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]','CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']"
```

This can become impractical and perhaps you prefer to provide an **input file** instead. Let's name this file `input.csv`.

{% code title="input.csv" %}
```bash
C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
```
{% endcode %}

The terminal command now becomes much cleaner:

```bash
# predict using an input file
ersilia api run -i input.csv
```

By default, predictions are returned in the standard **output** of the terminal. We favour the widely used **JSON format** because it offers great flexibility and interoperability. However, many of the model APIs return an output that can be naturally expressed in tabular format, for example, in a **CSV file**. If this is what you want, simply specify an output file with the `.csv` extension.

```bash
# save output in a CSV file
ersilia api run -i input.csv -o output.csv
```

{% hint style="info" %}
At the moment, the available formats are JSON (`.json`), CSV (`.csv`), TSV (`.tsv`) and HDF5 (`.h5`). The latter is appropriate for large-scale numerical data and is relevant for the lake of pre-computed predictions available in the [Isaura](https://github.com/ersilia-os/isaura) resource.
{% endhint %}

### Close model

Once you are done with predictions, it is advised to stop the model server:

```bash
# close model
ersilia close
```

### Delete model

If you are sure you don't want to use a model anymore, you may want to remove it from your computer. This includes deleting all model files and specific dependencies:

```bash
# delete model
ersilia delete retrosynthetic-accessibility
# or use the eos identifier
ersilia delete eos2r5a
```

## As a Python package

Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called `ErsiliaModel`:

```python
# import main class
from ersilia import ErsiliaModel
# instantiate the model
mdl = ErsiliaModel("retrosynthetic-accessibility")
```

Then, you can perform the same actions as in the CLI. To **serve**:

```python
# serve model
mdl.serve()
```

To make **predictions** for Halicin and Ibuprofen:

```python
# Halicin and Ibuprofen
input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
]
# predict
mdl.run(input)
```

To **close** the model:

```python
# close model
mdl.close()
```

### Using the `with` statement

A more concise way to run prediction would be to use the `with` clause:

```python
# use with statement
with ErsiliaModel("retrosynthetic-accessibility") as mdl:
    mdl.run(input)
```

## Using Ersilia through Colab

We have prepared a ready-to-go Google Colaboratory (Colab) notebook to run models and store predictions in Google Drive. Read more about Colab in our [training materials](../training-materials/google-colaboratory.md) and get started by clicking on the button Open in Colab from our [GitHub](https://github.com/ersilia-os/ersilia/blob/master/notebooks/ersilia-on-colab.ipynb).

{% hint style="warning" %}
Please note we do not extensively maintain the Colaboratory implementation and some models might not work. Please open an issue on GitHub if you encounter any problems.
{% endhint %}
