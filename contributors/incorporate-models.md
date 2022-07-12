---
description: This tutorial explains how to incorporate models in the Ersilia Model Hub
---

# Incorporate models

{% hint style="danger" %}
This page is **work in progress**!
{% endhint %}

## Anatomy of the Ersilia Model Template

Each model in the Ersilia Model Hub is contained within an individual GitHub repository. The [**Ersilia Model Template**](https://github.com/ersilia-os/eos-template) repository is stored as a GitHub Template, so you can [create a new repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template) based on it.

Below, we describe the template files in detail. Note that we only explain the files that you need to modify; other files, like `pack.py`, do not need any modification from the model contributor.

### The `eos` identifier

Each model in the Ersilia Model Hub has an Ersilia Open Source (EOS) identifier. This identifier determines the name of the GitHub repository containing the model:

```
https://github.com/ersilia-os/[EOS_IDENTIFIER]
```

The `eos` identifier follows this regular expression: `eos[1-9][a-z0-9]{3}`. That is:

* The `eos` prefix, plus...
* one digit  (`1-9`) (the `0` is reserved for test models), plus...
* three alphanumeric (`a-z` and `0-9`) characters.

{% hint style="success" %}
A list of `eos` identifiers needs is provided in the [Ersilia Model Hub Spreadsheet](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=drive\_web\&ouid=114775674178390159004). You can read about this Spreadsheet [here](model-selection.md).
{% endhint %}

### The [`README`](https://github.com/ersilia-os/eos-template/blob/main/README.md) file

The `README.md` file is where the basic information about the model is provided. It must include the following fields:

**Title:** a self-descriptive model title (less than 70 characters)

**Model Identifiers:** a set of codes that identify the model in the Ersilia CLI:

* Ersilia identifier (EOS ID): the `eos` identifier described above. Use the assigned identifier in the [Ersilia Model Hub Spreadsheet](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=sharing).
* Slug: a one-word or multi-word (linked by a hypen) human-readable identifier to be used as analternative to the EOS ID.
* Tags: labels to facilitate model search. For example, a model that predicts activity against malaria could have _P.Falciparum_ as tag. Select three relevant ones.

**Description**: minimum information about model type, results and the training dataset.

* Input: data format required by the model. Most chemistry related models, for example, will require molecules in SMILES format. If other input types, such as InChIKeys or protein sequences are accepted, specify them.
* Output: unit and description of the model result. It is essential to specify if the model gives back a probability or a specific measure. _For example, IC50_
* Model type: regression or classification
* Training set: number of compounds and link to the training dataset if available
* Mode of training: pretrained (the checkpoints where downloaded directly from a third party) retrained (the model was trained again using the same or a new dataset (please specify)) new model (if the model has been developed from scratch by Ersilia's contributors)

Results interpretation: provide a brief description of how to interpret the model results. _For example, in the case of a binary classification model for antimalarial activity based on experimental IC50, indicate the experimental settings (time of incubation, strain of parasite...) and the selected cut-off for the classification._

**Source Code:** this section must contain **all** relevant information about the model original authors, including a link to the publication if the model has been published in a peer reviewed journal or is in a preprint repositories, a link to the source code (typically, GitHub, GitLab or BitBucket) and a link to the model checkpoints directly, when available.

**License:** in addition to the `LICENSE` file, it is good practice to specify the Licenses in the README.md. All models in the Ersilia Model Hub are licensed under an open source license. Please make sure to abide by requirements of the original license when re-licensing or sub-licensing third-party author code. You can read more about how we deal with Open Source Licenses [here](https://ersilia.gitbook.io/ersilia-book/contributors/open-source-licences).

**History:** a short, numbered explanation of the modifications made in the source code, including the date of download and any steps taken to embed the model within the Ersilia Model Hub infrastructure.&#x20;

**About us:** all Ersilia repositories contain a final _About_ section, please keep the predefined version.

{% hint style="info" %}
The `README.md` file is written in [Markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) language. Please respect the hierarchy of headings provided in the template. Headings are specified with one or multiple `#` characters.
{% endhint %}

### The [`LICENSE`](https://github.com/ersilia-os/eos-template/blob/main/LICENSE) file

By default, all code written in contribution to Ersilia should be licensed under a [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). The main `LICENSE` file of the repository, therefore, will be a GPLv3 as specified by [GitHub.](https://docs.github.com/en/communities/setting-up-your-project-for-healthy-contributions/adding-a-license-to-a-repository)

However, the license notices for code developed by **third parties** must be kept in the respective folders where the third-party code is found. Include an explanation in the `README` file, for example:

> The GPLv3 license applies to all parts of the repository that are not externally maintained libraries. This repository uses the externally maintained library ChemProp library, located at `./model` and licensed under an [MIT License](https://github.com/ersilia-os/eos4e40/blob/main/model/LICENSE.md).

### The [`Dockerfile`](https://github.com/ersilia-os/eos-template/blob/main/Dockerfile) file

Ersilia uses a `Dockerfile` to specify installation instructions. The reason for this is that Docker provides the maximum level of isolation possible (i.e. a container), which may be needed to run models in some systems. However, in most practical scenarios, a Docker container will not be necessary and a Conda environment, or even a Virtualenv environment, will suffice. The Ersilia CLI will decide which isolation level to provide depending on the content of the `Dockerfile`:

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

The `Dockerfile` can contain as many `RUN` commands as necessary, between the `MAINTAINER` and the `WORKDIR` lines.

{% hint style="info" %}
The `Dockerfile` contains the installation instructions of the model. Therefore, the content of this file can be very variable, since each model will have its own dependencies.
{% endhint %}

### The [`model`](https://github.com/ersilia-os/eos-template/tree/main/model) folder

The `model` folder is the most important one. It contains two subfolders:

* `framework`: In this folder, we keep all the necessary code to run the AI/ML model (assuming dependencies are already installed).
* `checkpoints`: In this folder, we store the model data (pre-trained model parameters, scaling data, etc).

{% hint style="danger" %}
The `model` folder **should not** contain anything other than the `framework` and `checkpoints` subfolder. When the Ersilia CLI eventually fetches the model, it does a reorganization of the code and the only subfolders it keeps are these two. Any other file or folder at the `model/` directory level will be overlooked.
{% endhint %}

{% hint style="info" %}
Often, the separation between `framework` and `checkpoints` is not easy to determine. Sometimes, models provided from third-parties have model data embedded within the code or as part of the repository. In these cases, it is perfectly fine to keep model data in the `framework` subfolder, and leave the `checkpoints` subfolder empty.
{% endhint %}

The `framework` subfolder contains at least one Bash file, named `run_[API_NAME].sh`. Many models will have an API called `predict`, so `run_predict.sh`is frequently used. This file will run as follows:

```bash
bash run_predict.sh [DATA_FILE] [OUTPUT_FILE]
```

Unless strictly necessary, the `run_predict.sh` file should accept two and only two arguments, `DATA_FILE` and `OUTPUT_FILE`. In the current template, we provide the following example:

{% code title="run_predict.sh" %}
```bash
python src/main.py -i $1 -o $2
```
{% endcode %}

In this case, a Python file is executed, taking as input (`-i`) the `DATA_FILE` and giving as output (`-o`) the `OUTPUT_FILE`.

We now need to inspect the `main.py`file in more detail, since in this&#x20;

Unless strictly necessary, the `run_predict.sh` file should accept two and only two arguments, namely `DATA_FILE` and `OUTPUT_FILE`. In the current template we provide a simple example:

{% code title="bash_predict.sh" %}
```bash
FRAMEWORK_PATH = $(dirname -- "$(readlink -f "${BASH_SOURCE}")")

python {$FRAMEWORK_PATH}/src/main.py -i $1 -o $2
```
{% endcode %}

In this case, a Python script is run having as arguments an input (`-i`) and output (`-o`) files. Note that a `FRAMEWORK_PATH` variable is defined. This is used to store the directory name of the current script. The current template proposes the following script:

{% code title="src/main.py" %}
```python
import os
import csv
import os

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints directory
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# read checkpoints (here, simply an integer number)
ckpt = joblib.load(os.path.join(checkpoints_dir, "checkpoints.joblib"))

# model to be run (here, calculate the SMILES length and add ckpt to it)
def my_model(smiles_list, ckpt):
    return [len(smi)+ckpt for smi in smiles_list]
    
# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader) # skip header
    smiles_list = [r[0] for r in reader]
    
# run model
outputs = my_model(smiles_list, ckpt)

# write output in a .csv file
with open(output_file, "r") as f:
    writer = csv.writer(f)
    writer.writerow(["counts"]) # header
    for o in outputs:
        writer.writerow([o])
```
{% endcode %}

In this case, the model simply counts the number of characters in a SMILES string and adds a number to it.

The important steps of the script are:

1. Load model parameters.
2. Read input file.
3. Run predictions using the input file and the model parameters.or&#x20;
4. Write the output.

Most of the work of the model contributor will be working on this script. In the template, we provide a dummy model, which can be already defined within the script. Often, the model will be loaded from a third party Python library, or from a (cloned) repository placed in the same directory.

To summarize, in the template, we provide the a structure that follows this logic:

1. A `run_predict.sh` script executes a Python `main.py` script.
2. The `main.py` script:
   * Defines the model code.
   * Loads parameters from `checkpoints`.
   * Reads an input file containing SMILES (with header).
   * Runs a model that counts SMILES length and adds an integer defined by the parameters.
   * Writes an output file containing one column, i.e. the output value (with header)

{% hint style="info" %}
In the template, the example provided is very simple. Depending on the model being incorporated, the logic may be different. For example, many third party models already contain a command-line option, with a specific syntax. In these cases, you may want to write scripts to adapt the input and the output, and execute the model as-is.
{% endhint %}

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

This class is simply a wrapper for the AI/ML model. Typically, the central method of the `Model` class is the `predict` method.

```python
class Model (object):
    ...
    def predict(self, smiles_list):
        ...
```

In this case, the model takes as input a list of molecules represented as SMILES strings. This is the standard input type for [Type A](model-selection.md) models, focussed on chemistry data as input.

{% hint style="info" %}
You can always rename the `predict` method to something else if your model does not do predictions, strictly. For example, for some models it is more appropriate to rename to `calculate`.
{% endhint %}

{% hint style="success" %}
Multiple methods are allowed. For example, a model may have a `predict` and an `explain`method.
{% endhint %}

In its simplest form, the `Model` class just points Ersilia to the `model` directory and then creates a bash file to execute the necessary commands to run the model. So it is actually a very simple class, although it may look overwhelming at first. We break it down below:

First, a temporary directory is created:

```python
 class Model(object):
    ...
    def predict(self, smiles_list):
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
        data_file = os.path.join(tmp_folder, self.DATA_FILE)
        output_file = os.path.join(tmp_folder, self.OUTPUT_FILE)
        log_file = os.path.join(tmp_folder, self.LOG_FILE)
        ...
```

Then, a data file is created in the temporary directory. In this example, it is simply a one-column `csv` file having a header (`smiles`) and a list of molecules in SMILES format (one per row):

```python
class Model(object):
    ...
    def predict(self, smiles_list):
        ...
        with open(data_file, "w") as f:
            f.write("smiles"+os.linesep)
            for smiles in smiles_list:
                f.write(smiles+os.linesep)
        ...
```

Now we already have the input file ready to submit to the corresponding `run_predict.sh`file in the `model/framework/` directory, as specified above. The following creates a dummy bash script in the temporary directory and runs the command from there. The output is saved in the temporary directory too.

```python
class Model(object):
    ...
    def predict(self, smiles_list):
        ...
        run_file = os.path.join(tmp_folder, self.RUN_FILE)
        with open(run_file, "w") as f:
            lines = [
                "bash {0}/run_predict.sh {1} {2}".format( # <-- EDIT: match method name (run_predict.sh, run_calculate.sh, etc.)
                    self.framework_dir,
                    data_file,
                    output_file
                )
            ]
            f.write(os.linesep.join(lines))
        cmd = "bash {0}".format(run_file)
        with open(log_file, "w") as fp:
            subprocess.Popen(
                cmd, stdout=fp, stderr=fp, shell=True, env=os.environ
            ).wait()
        ...
```

The last step is to read from the output in the temporary directory and return it in a JSON-serializable format. The output in the example is a `csv` table, with one or multiple columns, containing numeric data. The table has a header, which is read and saved as metadata.

```python
class Model(object):
    ...
    def predict(self, smiles_list):
        ...
        with open(output_file, "r") as f:
            reader = csv.reader(f)
            h = next(reader)
            R = []
            for r in reader:
                R += [{"outcome": [Float(x) for x in r]}]
        meta = {
            "outcome": h
        }
        result = {
            "result": R,
            "meta": meta
        }
        return result
```

You will see that, in the template, pointers to potential edits are highlighted with the tag `# EDIT` . Necessary edits relate to the appropriate naming of the API, the format of the input data, or the serialization to JSON format from the output data.

{% hint style="info" %}
Advanced contributors may want to modify the `Model` class to load a model in-place (for example, a Scikit-Learn model) instead of executing a bash command in the `model/framework/` directory.
{% endhint %}

#### The `Artifact` class

This class mirrors BentoML artifacts. It simply contains `load`, `save`, `get and pack` functionalities:

```python
class Artifact(BentoServiceArtifact):
    ...
    def pack(self, model):
        ...
    def load(self, path):
        ...
    def get(self):
        ...
    def save(self, dst):
        ...
```

You **don't have to modify** this class.

#### The `Service` class

This class is used to create the service. The service at least one API, typically `predict`:&#x20;

```python
@artifacts([Artifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def predict(self, input: List[JsonSerializable]):
        input = input[0]
        smiles_list = [inp["input"] for inp in input]
        output = self.artifacts.model.predict(smiles_list)
        return [output]
```

The `Service`class can have multiple APIs, each of them specified with the `@api` decorator. By default, Ersilia works with JSON inputs, which are deserialized as a SMILES list inside the API.

You can **rename the API** (for example, to `calculate`), following the `# EDIT` tags in the template.

## Steps for model incorporation

In this tutorial, we will follow the example of a very simple model, related synthetic accessibility scoring.

#### Find model code and parameters

XX

#### Run the code _outside_ Ersilia

XX

## TL;DR

In summary, the steps to incorporate a model to the Ersilia Model Hub are the following. We are assuming that the model can be installed in a Conda environment.

1. Download the model from a third party repository to your local machine.
2. Install the model in a dedicated Conda environment and make sure you can run it.
3. Create a new GitHub repository from the `eos-template`. Name this repository with the EOS ID.
4. Clone the new repository.
5. Place model code in `model/framework` and model parameters in `model/checkpoints`.
6. Write the necessary code to provide a `run_predict.sh` that simply takes one input file and produces one output file. Be sure to use absolute paths throughout.
7. Edit the `Dockerfile` file to reflect the installation steps followed in 2.
8. Edit the `service.py` file as needed.
9. Make sure that `.gitattributes` specifies the format of your model parameters.
10. Push changes to the model repository.
11. Activate the Ersilia CLI and fetch the model (in verbose mode).
12. Serve and run the model.



