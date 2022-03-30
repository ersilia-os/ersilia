# Antibiotic activity prediction

In this case example, we show how to run predictions based on the AI/ML model developed in the paper _A Deep Learning Approach to Antibiotic Drug Discovery_, by [Stokes et al. Cell (2020)](https://www.cell.com/cell/fulltext/S0092-8674\(20\)30102-1?\_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420301021%3Fshowall%3Dtrue). In this study, a deep learning model was trained using the excellent [Chemprop](https://github.com/chemprop/chemprop) tool. Halicin, a drug originally researched for the treatment of diabetes, was predicted to have broad antibacterial activity.

## Command-line interface

We provide a command-line interface (CLI) to interact with the Ersilia Model Hub. To check the available commands, simply type:

```bash
# list available commands
ersilia --help
```

### Browse the model catalog

You can explore our catalog of models. The following will return a list of models currently available in our remote repositories.

```bash
# catalog of models
ersilia catalog
```

{% hint style="info" %}
The Ersilia Model Hub is growing continuously to fulfil the needs of the community. Please do not hesitate to request new models! Just reach out to us, we will be happy to assist: [hello@ersilia.io](mailto:hello@ersilia.io)
{% endhint %}

In the Ersilia Model Hub, the antibiotic activity prediction model has the **identifier** `eos4e40`  and the **slug** `chemprop-antibiotic`. You can use either one to refer to this model all of the commands below. For example, you can get more information through the model card:

```bash
# display model card using slug...
ersilia card chemprop-antibiotic
# ... or using ersilia identifier
ersilia card eos4e40
```

### Fetch model and install it locally

The first step is to download the model to your local device and install it along with its dependencies. By default, a `~/eos` directory (for Ersilia Open Source) will be created in your `HOME`. This folder will contain all fetched models along with additional files to manage the AI/ML content available locally.

To download and install the antibiotic activity prediction model, simply use the `fetch` command:

```bash
# fetch model from remote repository
ersilia fetch chemprop-antibiotic
```

{% hint style="warning" %}
This model weighs 850 MB, so fetching will take a while. In order to provide robust predictions, the original paper released an ensemble of 20 individual models. We do not want to compromise the accuracy of the model, this is why the full ensemble is downloaded. Our team is currently developing a **lightning** library that will provide lighter versions of the models for easier testing and deployment.
{% endhint %}

{% hint style="info" %}
We do our best to keep the user away from the [dependency hell](https://en.wikipedia.org/wiki/Dependency\_hell). Models are **automatically** installed with the necessary degree of isolation from the system. While some models have no dependencies at all and can be run using the system Python installation, others need to be containerized using Docker.&#x20;
{% endhint %}

### Serve model

Once the model has been fetched, it should be ready to be used. A model in the Ersilia Model Hub can be thought of as a set of APIs. You can serve the model like this:

```bash
# serve model
ersilia serve chemprop-antibiotic
```

A URL will be prompted as well as a process id (PID). These can be relevant if you are an advanced user and want to have low-level control of the tool. The most important is, however, the list of **available APIs**. In this case, we want to infer antibiotic activity through the `predict` API.

{% hint style="info" %}
The `predict` API is obviously one of the most ubiquitous throughout our catalog of models. Other common APIs are `transform` and `interpret`.
{% endhint %}

### Make predictions

The antibiotic activity prediction model takes **chemical structures** as input and provides an activity score (ranging from 0 to 1) under a cutoff of _E.coli_ growth inhibition of 50 uM.

Ideally, in the chemistry models, the input molecules are specified as **SMILES** strings. SMILES strings can be easily found online. For instance, we can find Halicin in [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/Halicin#section=Canonical-SMILES) and then predict its antimicrobial activity as follows:

```bash
# Halicin
ersilia api predict -i "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]"
```

{% hint style="warning" %}
It is also possible to use [InChIKey](https://pubchem.ncbi.nlm.nih.gov/compound/Halicin#section=InChI-Key) or even molecule name (through the [Chemical Identifier Resolver](https://cactus.nci.nih.gov/chemical/structure)) instead of SMILES. Ersilia will take care of this automatically. However, please take into account that this requires an internet connection and will slow down the process, as requests to external tools are necessary.
{% endhint %}

You can make **multiple predictions** in batch mode. This is typically much faster than running predictions one by one in a loop. For instance, we can predict the antimicrobial activity of Halicin and [Ibuprofen](https://pubchem.ncbi.nlm.nih.gov/compound/Ibuprofen#section=Canonical-SMILES). We don't expect Ibuprofen to be active.

```bash
# Halicin and Ibuprofen
ersilia api predict -i "['C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]','CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']"
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
ersilia api predict -i input.csv
```

By default, predictions are returned in the standard **output** of the terminal. We favour the widely used **JSON format** because it offers great flexibility and interoperability. However, many of the model APIs return an output that can be naturally expressed in tabular format, for example, in a **CSV file**. If this is what you want, simply specify an output file with the `.csv` extension.

```bash
# save output in a CSV file
ersilia api predict -i input.csv -o output.csv
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
ersilia delete chemprop-antibiotic
```

## As a Python package

Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called `ErsiliaModel`:

```python
# import main class
from ersilia import ErsiliaModel
# instantiate the model
mdl = ErsiliaModel("chemprop-antibiotic")
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
mdl.predict(input)
```

To **close** the model:

```python
# close model
mdl.close()
```

### Using the 'with' statement

A more concise way to run prediction would be to use the `with` clause:

```python
# use with statement
with ErsiliaModel("chemprop-antibiotic") as mdl:
    mdl.predict(input)
```
