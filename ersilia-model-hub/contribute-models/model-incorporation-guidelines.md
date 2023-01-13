---
description: This tutorial explains how to incorporate models in the Ersilia Model Hub
---

# Model Incorporation Guidelines

{% hint style="danger" %}
We are currently working on a model incorporation pipeline. Content in this page is therefore slightly outdated.
{% endhint %}

## Anatomy of the Ersilia Model Template

Each model in the Ersilia Model Hub is contained within an individual GitHub repository. The [**Ersilia Model Template**](https://github.com/ersilia-os/eos-template) repository is stored as a GitHub Template, so you can [create a new repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template) based on it.

Below, we describe the template files in detail. Note that we only explain the files that you need to modify; other files, like `pack.py`, do not need modification from the model contributor.

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
The list of `eos` identifiers is pre-coded, please reach out to one of the Ersilia Model Hub managers, and do not come up with your own identifier
{% endhint %}

### The [`README`](https://github.com/ersilia-os/eos-template/blob/main/README.md) file

The `README.md` file is where we give basic information about the model. It must include the following fields:

**Title:** a self-descriptive model title (less than 70 characters)

**Model Identifiers:** a set of codes that identify the model

* Ersilia identifier (EOS ID): the `eos` identifier described above. Use the assigned identifier.
* Slug: a one-word or multi-word (linked by a hypen) human-readable identifier to be used as an alternative to the EOS ID.
* Tags: labels to facilitate model search. For example, a model that predicts activity against malaria could have _Plasmodium falciparum_ as tag. Select three relevant tags.

**Description**: minimum information about model type, results and the training dataset.

* Input: data format required by the model. Most chemistry related models, for example, will require molecules in SMILES format. If special input types are required, please specify them.
* Output: description of the model result. It is important to be precise in this description. Is the model providing a probability? Is it a score? Is it single-output or multi-output? etc.
* Model type: regression, classification, embedding...
* Training set: number of compounds and link to the training dataset, if available
* Mode of training: pretrained (the checkpoints where downloaded directly from a third party) retrained (the model was trained again using the same or a new dataset), new (if the model has been developed from scratch by Ersilia's contributors)

{% hint style="info" %}
Some contributors may find it difficult to come up with a good description for the model. You can find some inspiration in [Semantic Scholar](https://semanticscholar.org). This portal provides an AI-based **TL;DR** short description of many indexed papers. \
\
You can also try [ChatGPT](https://openai.com/blog/chatgpt/)!
{% endhint %}

**Results interpretation:** provide a brief description of how to interpret the model results. For example, in the case of a binary classification model for antimalarial activity based on experimental IC50, indicate the experimental settings (time of incubation, strain of parasite...) and the selected cut-off for the classification.

**Source code:** this section must contain _all_ relevant information about the original authors of the model, including a link to the publication if the model has been published in a peer reviewed journal or is in a preprint repository, a link to the source code (typically, GitHub, GitLab or BitBucket) and a link to the model checkpoints directly, when available.

**License:** in addition to the `LICENSE` file, it is good practice to specify the licenses in the `README` file. All models in the Ersilia Model Hub are licensed under an open source license. Please make sure to abide by requirements of the original license when re-licensing or sub-licensing third-party author code. You can read more about how we deal with Open Source Licenses [here](https://ersilia.gitbook.io/ersilia-book/contributors/open-source-licences).

**History:** a short, numbered explanation of the modifications made in the source code, including the date of download and any steps taken to embed the model within the Ersilia Model Hub infrastructure.&#x20;

**About us:** all Ersilia repositories contain a final _About_ section, please keep the predefined version.

{% hint style="info" %}
The `README` file is written in [Markdown](https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax) language. Please respect the hierarchy of headings provided in the template. Headings are specified with one or multiple `#` characters.
{% endhint %}

### The [`LICENSE`](https://github.com/ersilia-os/eos-template/blob/main/LICENSE) file

By default, all code written in contribution to Ersilia should be licensed under a [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). The main `LICENSE` file of the repository, therefore, will be a GPLv3 as specified by [GitHub.](https://docs.github.com/en/communities/setting-up-your-project-for-healthy-contributions/adding-a-license-to-a-repository)

However, the license notices for code developed by **third parties** must be kept in the respective folders where the third-party code is found. Include an explanation in the `README` file, for example:

> The GPLv3 license applies to all parts of the repository that are not externally maintained libraries. This repository uses the externally maintained library ChemProp library, located at `./model` and licensed under an [MIT License](https://github.com/ersilia-os/eos4e40/blob/main/model/LICENSE.md).

### The [`Dockerfile`](https://github.com/ersilia-os/eos-template/blob/main/Dockerfile) file

Ersilia uses a `Dockerfile` file to specify installation instructions. The reason for this is that Docker provides the maximum level of isolation possible (i.e. a container), which may be needed to run models in some systems. However, in most practical scenarios, a Docker container will not be necessary and a Conda environment, or even a Virtualenv environment, will suffice. The Ersilia CLI will decide which isolation level to provide depending on the content of the `Dockerfile`:

* If only `pip install` commands are specified, Virtualenv will be used.
* If only `pip install` and `conda install` commands are specified, Conda will be used.
* If other commands are specified (e.g. `sudo apt-get`), Docker will be used.

The `Dockerfile` available in the Ersilia Model Template is the following:

{% code title="Dockerfile" %}
```docker
FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit=2020.03
RUN pip install joblib==1.1.0

WORKDIR /repo
COPY . /repo
```
{% endcode %}

In this case, a Conda environment will be preferentially used to isolate the model. The first line of the `Dockerfile` indicates that this Conda environment will have **BentoML 0.11.0** installed on **Python 3.7**.

In this example, the `rdkit` library will be installed using `conda`, and `joblib` will be installed using `pip`.

The `Dockerfile` can contain as many `RUN` commands as necessary, between the `MAINTAINER` and the `WORKDIR` lines.

{% hint style="warning" %}
The `Dockerfile` contains the installation instructions of the model. Therefore, the content of this file can be very variable, since each model will have its own dependencies.
{% endhint %}

### The [`model`](https://github.com/ersilia-os/eos-template/tree/main/model) folder

The `model` folder is the most important one. It contains two subfolders:

* `framework`: In this folder, we keep all the necessary code to run the model (assuming dependencies are already installed).
* `checkpoints`: In this folder, we store the model data (pretrained model parameters, scaling data, etc).

{% hint style="danger" %}
The `model` folder **should not** contain anything other than the `framework` and `checkpoints` subfolder. When the Ersilia CLI eventually fetches the model, it does a reorganization of the code and the only subfolders it keeps are these two. Any other file or folder at the `model/` directory level will be overlooked.
{% endhint %}

{% hint style="info" %}
Often, the separation between `framework` and `checkpoints` is not easy to determine. Sometimes, models obtained from third parties have model data embedded within the code or as part of the repository. In these cases, it is perfectly fine to keep model data in the `framework` subfolder, and leave the `checkpoints` subfolder empty.
{% endhint %}

The `framework` subfolder contains at least one Bash file, named `run_[API_NAME].sh`. Many models will have an API called `predict`, so `run_predict.sh`is frequently used. This file will run as follows:

```bash
bash run_predict.sh [FRAMEWORK_DIR] [DATA_FILE] [OUTPUT_FILE]
```

Unless strictly necessary, the `run_predict.sh` file should accept three and only three arguments, namely `FRAMEWORK_DIR`, `DATA_FILE` and `OUTPUT_FILE`. In the current template, we provide the following example:

{% code title="run_predict.sh" %}
```bash
python $1/code/step.py -i $2 -o $3
```
{% endcode %}

In this case, a Python file located in the `[FRAMEWORK_DIR]/src` folder is executed, taking as input (`-i`) the `DATA_FILE` and giving as output (`-o`) the `OUTPUT_FILE`.

To understand this further, we now need to inspect the step`.py`file in more detail. The current template proposes the following script:&#x20;

{% code title="code/step.py" %}
```python
# imports
import os
import csv
import joblib
import sys
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints directory
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# read checkpoints (here, simply an integer number: 42)
ckpt = joblib.load(os.path.join(checkpoints_dir, "checkpoints.joblib"))

# model to be run (here, calculate the Molecular Weight and add ckpt (42) to it)
def my_model(smiles_list, ckpt):
    return [MolWt(Chem.MolFromSmiles(smi))+ckpt for smi in smiles_list]
    
# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader) # skip header
    smiles_list = [r[0] for r in reader]
    
# run model
outputs = my_model(smiles_list, ckpt)

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["value"]) # header
    for o in outputs:
        writer.writerow([o])
```
{% endcode %}

In this case, the model simply calculates the molecular weight and adds a number to it.

The important steps of the script are:

1. Load model parameters.
2. Read input file.
3. Run predictions using the input file and the model parameters.
4. Write the output.

Most of the work of the model contributor will be to work on this or similar scripts. In the template, we provide a dummy model (i.e. add a fixed value to the molecular weight). This dummy model can can be already defined within the script (`my_model`). However, in real world cases, the model will most likely be loaded from a third party Python library, or from a (cloned) repository placed in the same directory.

To summarize, in the template, we provide a structure that follows this logic:

1. The `run_predict.sh` script executes the Python `main.py` script.
2. The `step.py` script:
   * Defines the model code.
   * Loads parameters from `checkpoints`.
   * Reads an input file containing SMILES (with header).
   * Runs a model that calculates molecular weight and adds an integer defined by the parameters.
   * Writes an output file containing one column corresponding to the output value (with a header).

{% hint style="info" %}
In the template, the example provided is very simple. Depending on the model being incorporated, the logic may be different. For example, many third party models already contain a command-line option, with a specific syntax. In these cases, you may want to write scripts to adapt the input and the output, and then execute the model as-is.

Each script will be one step.py file, we can create as many as necessary and rename them appropriately (see below for examples)
{% endhint %}

### The [`.gitattributes`](https://github.com/ersilia-os/eos-template/blob/main/.gitattributes) file

We use Git LFS to store large files (over 100 MB). Typically, these files are model parameters. Files to be stored in Git LFS should be specified in the `.gitattributes` file. The current file will store in Git LFS all files in `csv`, `h5`, `joblib`, `pkl`, `pt` and `tsv` format.

```
*.csv filter=lfs diff=lfs merge=lfs -text
*.h5 filter=lfs diff=lfs merge=lfs -text
*.joblib filter=lfs diff=lfs merge=lfs -text
*.pkl filter=lfs diff=lfs merge=lfs -text
*.pt filter=lfs diff=lfs merge=lfs -text
*.tsv filter=lfs diff=lfs merge=lfs -text
```

You have to edit the `.gitattributes` file to ensure that all large files in your model are stored in Git LFS. The `git lfs track` command automatically updates this file.

### The [`service`](https://github.com/ersilia-os/eos-template/blob/main/src/service.py) file

The service file is located in `src/service.py`. It contains the necessary code to facilitate model bundling with BentoML.

There are three main classes in the `service` file, namely `Model`, `Artifact` and `Service`.

#### The `Model` class

This class is simply a wrapper for the AI/ML model. Typically, when incorporating **external** (type 1) models, the `run_predict.sh` script will already capture the logic within the `Model` class, in which case the `Model` class is simply redundant. However, when incorporating **internally developed** (types 2 and 3) models into the hub, we can make use of the artifacts for standard modeling frameworks (e.g. sklearn, PyTorch, and Keras) provided by BentoML, and the `Model` class becomes necessary for BentoML compatibility. Hence, the `Model` class enables generalization between these types of model additions.

Typically, the central method of the `Model` class is the `predict` method.

```python
class Model (object):
    ...
    def predict(self, smiles_list):
        ...
```

In this case, the model takes as input a list of molecules represented as SMILES strings. This is the standard input type for [Type A](broken-reference) models, focussed on chemistry data as input.

{% hint style="info" %}
You can always rename the `predict` method to something else if your model does not do predictions, strictly. For example, for some models it is more appropriate to rename this method to `calculate`.
{% endhint %}

{% hint style="success" %}
Multiple methods are allowed. For example, a model may have a `predict` and an `explain`method.
{% endhint %}

In its simplest form, the `Model` class just points Ersilia to the `model` directory and then creates a Bash file to execute the necessary commands to run the model. It is actually a very simple class, although it may look overwhelming at first. We break it down below:

First, a temporary directory is created:

```python
 class Model(object):
    ...
    def predict(self, smiles_list):
        tmp_folder = tempfile.mkdtemp(prefix="eos-")
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

Now we already have the input file of the `run_predict.sh`script, located in the `model/framework/` directory, as specified [above](model-incorporation-guidelines.md#the-model-folder). The following creates a dummy Bash script in the temporary directory and runs the command from there. The output is saved in the temporary directory too. [Remember](model-incorporation-guidelines.md#the-model-folder) that the `run_predict.sh` script expects three arguments, `FRAMEWORK_DIR`_,_ `DATA_FILE` and `OUTPUT_FILE`.

```python
class Model(object):
    ...
    def predict(self, smiles_list):
        ...
        run_file = os.path.join(tmp_folder, self.RUN_FILE)
        with open(run_file, "w") as f:
            lines = [
                "bash {0}/run_predict.sh {0} {1} {2}".format(
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
        shutil.rmtree(tmp_folder)
        return result
```

You will see that, in the template, pointers to potential edits are highlighted with the tag `# EDIT` . Necessary edits relate to the appropriate naming of the API, the format of the input data, or the serialization to JSON format from the output data.

{% hint style="info" %}
Advanced contributors may want to modify the `Model` class to load a model in-place (for example, a Scikit-Learn model) instead of executing a Bash command in the `model/framework/` directory.
{% endhint %}

#### The `Artifact` class

This class mirrors [BentoML artifacts](https://docs.bentoml.org/en/0.13-lts/api/index.html). It simply contains `load`, `save`, `get and pack` functionalities:

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

This class is used to create the service. The service exposes at least one API, typically `predict`:&#x20;

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

The `Service`class can have multiple APIs, each of them specified with the `@api` decorator. By default, Ersilia works with JSON inputs, which are deserialized as a SMILES list inside the API, in this case.

You can **rename the API** (for example, to `calculate`), following the `# EDIT` tags in the template.

## Steps for model incorporation

Now that we have an idea of the contents of the [Ersilia Model Template](https://github.com/ersilia-os/eos-template), we will follow the example of a simple but widely used model to calculate the **synthetic accessibility** of small molecule compounds. Synthetic accessibility measures the feasibility of synthesizing a molecule in the laboratory. In [2009, Peter Ertl presented the synthetic accessiblity (SA) score](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8), based on measures of molecular complexity and occurrence of certain fragments in the small molecule structure. High (greater than 6) SA scores denote difficult-to-synthesize molecules, and low (lower than 3) SA scores suggest that the molecule will be easy to synthesize.

### Include the model to the Ersilia Model Hub AirTable

The [Ersilia CLI](https://github.com/ersilia-os/ersilia) accesses data contained in the [Ersilia Model Hub AirTable](https://airtable.com/shrUcrUnd7jB9ChZV/tblZGe2a2XeBxrEHP) database. Thus, the first step is to include the model entry in the database. You can follow the instruction in the [model selection](broken-reference) page. You will see that, in the Ersilia Model Hub, the current model has the EOS identifier `eos-9ei3` and the slug `sa-score`.

{% hint style="info" %}
Please contact **@Miquel** if your model is not included in the AirTable database.
{% endhint %}

#### Read the publication

It is important that you read the original publication in order to understand the training data, the limitations of the model, the validation performed and the characteristic of the algorithm, among other details. In [this case](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8), this is a classic (old) publication from a Novartis team. They analyzed small fragments in the PubChem database and devised a molecular complexity score that takes into account molecule size, presence of non-standard structural features, presence of large rings, etc. The authors validated the model by comparing their SA score with synthetic feasibility as estimated by medicinal chemists for a set of 40 molecules.

#### Find model code and parameters

Code to calculate the SA score does not seem to be available from the publication. Fortunately, though, the RDKit library, in its contributions module, contains an implementation of the SA score. The code can be found [here](https://github.com/rdkit/rdkit/tree/master/Contrib/SA\_Score). This RDKit-based implementation was developed in 2013 by Peter Ertl and Greg Laundrum.

{% hint style="success" %}
Both the link to the code and to the original publication are accessible from the Ersilia Model Hub AirTable database.
{% endhint %}

### Run the code outside Ersilia

Before incorporating the `sa-score` model to the Ersilia Model Hub, we need to make sure that we can actually run the code provided by the third party. [In this case](https://github.com/rdkit/rdkit/tree/master/Contrib/SA\_Score), upon quick inspection, two elements seem to be central in the repository:

* The `sascorer.py` script, containing the main code. We can consider this file to be the **model code**.
* The `fpscores.pkl.gz` compressed file, containing pre-calculated fragment scores. In this simple case, we can consider this file to be the **model parameters**.

#### Create a conda environment

No installation instructions are provided for this model. However, the `sascorer.py` file `import` statements indicate that, at least, `rdkit` is necessary. We can create a Conda environment and install the `rdkit` as follows:

```bash
conda create -n sa-score python=3.7
conda activate sa-score
conda install -c conda-forge rdkit=2021.03
```

#### Download code and parameters

We can download code and parameters of the directly from the RDKit contributions repository. Here, we store these data in a folder named `SA_Score` located in the `~/Desktop`.

```bash
cd ~/Desktop
mkdir SA_Score
cd SA_Score
wget https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score/sascorer.py
wget https://raw.githubusercontent.com/rdkit/rdkit/master/Contrib/SA_Score/fpscores.pkl.gz
```

{% hint style="info" %}
Often, the model will be available as a full GitHub repository. In these cases, you can simply clone the repository.
{% endhint %}

#### Test the model

Inspection of the `sascorer.py` file indicates that we can run this script from the terminal:

{% code title="sascorer.py" %}
```python
...
if __name__ == '__main__':
    import sys
    import time

    t1 = time.time()
    readFragmentScores("fpscores")
    t2 = time.time()

    suppl = Chem.SmilesMolSupplier(sys.argv[1])
    t3 = time.time()
    processMols(suppl)
    t4 = time.time()

    print('Reading took %.2f seconds. Calculating took %.2f seconds' % ((t2 - t1), (t4 - t3)),
          file=sys.stderr)
```
{% endcode %}

The `Chem.SmilesMolSupplier` takes a file containing SMILES strings and identifiers, separated by a space character. The file expects a header. Let's create a file with a few molecules. We looked for examples in [DrugBank](https://go.drugbank.com/):

{% code title="molecules.smi" %}
```
smiles identifier
CC(=O)NC1=CC=C(O)C=C1 mol1
CN(C)CCC1=CNC2=C1C=C(CS(=O)(=O)N1CCCC1)C=C2 mol2
CCN(CC)CC1=C(O)C=CC(NC2=C3C=CC(Cl)=CC3=NC=C2)=C1 mol3
```
{% endcode %}

Now let's see if the model works as expected:

```bash
python sascorer.py molecules.smi
```

It does work! The output in the terminal looks like this:

```
smiles	Name	sa_score
CC(=O)Nc1ccc(O)cc1	mol1	1.407299
CN(C)CCc1c[nH]c2ccc(CS(=O)(=O)N3CCCC3)cc12	mol2	2.319770
CCN(CC)Cc1cc(Nc2ccnc3cc(Cl)ccc23)ccc1O	mol3	2.249844
Reading took 0.21 seconds. Calculating took 0.00 seconds
```

{% hint style="info" %}
Many repositories give a clear description of the expected input format. For the SA scorer, the expected input was not clearly specified, and previous knowledge of the `Chem.SmilesMolSupplier` method was necessary.
{% endhint %}

### Create the model repository from the Ersilia Model Template

Now that we know that the code can run in our local machine, we can create the corresponding respository in the Ersilia GitHub organization.

#### Use the eos-template

The easiest way to create the repository is to visit the `eos-template` [GitHub repository page](https://github.com/ersilia-os/eos-template). In green (top-right), you will find the **Use this template** button.

A page to create a new repository from `eos-template` will open:

* Set the **Owner** to `ersilia-os`.
* In the **Repository name**, enter the EOS identifier. In this case, `eos9ei3`.
* Use the **Description** provided in the [Ersilia Model Hub AirTable](https://airtable.com/shrUcrUnd7jB9ChZV/tblZGe2a2XeBxrEHP).
* Keep the repository **Public**.

Click **Create repository from template**.

#### Clone the new repository

We can now clone the repository in our local machine. Let's do it in the `~/Desktop`:

```bash
cd ~/Desktop
git clone https://github.com/ersilia-os/eos9ei3.git
cd eos9ei3
```

### Migrate code and parameters

Let's now place the code and the parameters in the `model` folder (in the `framework` and `checkpoints` subfolders, respectively):

```bash
cd ~/Desktop
cp SA_Score/sascorer.py eos9ei3/model/framework/code/.
cp SA_Score/fpscores.pkl.gz eos9ei3/model/checkpoints/.
```

{% hint style="info" %}
The checkpoints contains a dummy checkpoints.joblib that can be deleted
{% endhint %}

{% hint style="danger" %}
Note that here we are migrating code and parameters to different folders. This may cause critical errors if code expects to find parameters at a certain relative location. Try to locate the pointers to the model parameters and change the paths. Only, and exceptionally, if the model architecture is too complex we can keep code and parameters in the /framework folder. Please ask for permission to Ersilia's team before doing it.
{% endhint %}

### Write framework code

Now it is time to write some code. Here we will follow the description of the `model` folder [given above](model-incorporation-guidelines.md#the-model-folder).

#### Write input and output adapters

The `eos-template` provides an exemplary `step.py` that is not useful here. We will actually need to use three steps:

1. Input adapter
2. SAScorer (in this case, the `sascorer.py` which is already written)
3. Output adapter

```bash
cd model
rm framework/code/step.py
```

By default, for chemical compound inputs, Ersilia uses single-column files with a header (see the `service.py` file [above](model-incorporation-guidelines.md#the-service-file)). However, the `sascorer.py` [expects](model-incorporation-guidelines.md#test-the-model) a two-column file. Let's write an **input** adapter:

{% code title="code/input_adapter.py" %}
```python
import sys
import csv

input_file = sys.argv[1]
smiles_list = []
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader) # skip header
    for r in reader:
        smiles_list += [r[0]]

with open("tmp_input.smi", "w") as f:
    writer = csv.writer(f, delimiter=" ")
    writer.writerow(["smiles", "identifier"]) # header
    for i, smi in enumerate(smiles_list):
        writer.writerow([smi, "mol{0}".format(i)])
```
{% endcode %}

The script creates an intermediate `tmp_input.smi` file that can be used as input for `sascorer.py`.

Likewise, the **output** is expected to be just one column containing the SA score, in this case. The output provided by `sascorer.py` has three columns (tab-separated), so we need to adapt it.

{% code title="code/output_adapter.py" %}
```python
import sys
import csv

output_file = sys.argv[1]

sascores = []
with open("tmp_output.csv", "r") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader) # skip header
    for r in reader:
        sascores += [r[-1]]

with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["sa_score"]) # header
    for sa in sascores:
        writer.writerow([sa])
```
{% endcode %}

Note that we are reading from a `tmp_output.csv`. We then write a one-column output.

#### Make sure that parameters are read

So far, we haven't pointed to the model parameters. When [migrating code and parameters](model-incorporation-guidelines.md#migrate-code-and-parameters), we separated the `sascore.py` file and the `fpscores.pkl.gz` file.

Let's inspect `sascore.py` to understand how parameters are read. There is a `readFragmentScore` function that does this job. We need to modify it to point to the `checkpoints` folder:

{% code title="code/sascorer.py" %}
```python
...
def readFragmentScores(name='fpscores'):
    import gzip
    global _fscores
    # generate the full path filename:
    if name == "fpscores":
        name = op.join(op.join(op.dirname(__file__), "..", "..", "checkpoints"), name)
    data = pickle.load(gzip.open('%s.pkl.gz' % name))
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict
...
```
{% endcode %}

#### Write the `run_predict.sh` file

We now have the input adapter, the model code and parameters, and the output adapter. Let's simply write this pipeline in the `run_predict.sh` file:

{% code title="run_predict.sh" %}
```bash
python $1/code/input_adapter.py $2
python $1/code/sascorer.py tmp_input.smi > tmp_output.csv
python $1/code/output_adapter.py $3
rm tmp_input.smi tmp_output.csv
```
{% endcode %}

### Run predictions

Let's now check that the scripts run as expected. Eventually, Ersilia will run this code from an arbitrary location, so it is best to test it outside the `framework` folder. We can create an input file in the `~/Desktop`.

{% code title="molecules.csv" %}
```
smiles
[H][C@]12SC(C)(C)[C@@H](N1C(=O)[C@H]2NC(=O)[C@H](N)C1=CC=CC=C1)C(O)=O
NC1=CC(=CNC1=O)C1=CC=NC=C1
ClC1=CC=C2N=C3NC(=O)CN3CC2=C1Cl
COC1=CC=C(C=C1)C1=CC(=S)SS1
COC1=CC=C(C=C1)C(=O)CC(=O)C1=CC=C(C=C1)C(C)(C)C
```
{% endcode %}

To test the model, we simply have to execute the `run_predict.sh` script.

```
cd ~/Desktop
FRAMEWORK_PATH="eos9ei3/model/framework/"
bash $FRAMEWORK_PATH/run_predict.sh $FRAMEWORK_PATH molecules.csv output.csv
```

You should get the `output.csv` file in your `~/Desktop`. The output contains five predictions, corresponding to the five molecules in the `molecules.csv` file.

{% code title="output.csv" %}
```
sa_score
3.527006
2.344611
2.957539
2.535359
1.798579
```
{% endcode %}

### Edit the `service.py` file, if necessary

The `service.py` file provided by default in the template manages chemistry inputs and expects tabular outputs. Therefore, in principle, you do not have to modify this file.

{% hint style="info" %}
Modifying the `service.py` file is intended for advanced users only. Please use the Slack `#internships` channel if you think your model of interest requires modification of this file.
{% endhint %}

### Edit the `Dockerfile` file

The `Dockerfile` file should include all the installation steps that you run after creating the working Conda environment. In the case of `sa-score`, we only installed RDKit:

{% code title="Dockerfile" %}
```docker
FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit=2021.03.4

WORKDIR /repo
COPY . /repo
```
{% endcode %}

### Write the `README` file

Don't forget document the model. Read the [instructions to write the `README` file](model-incorporation-guidelines.md#the-readme-file) page. Feel free to ask for help in the Slack `#internships` channel.

### Commit changes to the repository

We are ready to commit changes and push them to the Ersilia Model Hub.

#### Check the `.gitattributes` file

If you have large files, you will have to track them with Git LFS. The template already provides a collection of common extensions to track with Git LFS (see `.gitattributes` file). However, it is possible that your parameters do not have any of these extensions.

Let's track our parameters, even if the file is not particularly big:

```bash
cd ~/Desktop/eos9ei3
git lfs track "*.pkl.gz"
git add .gitattributes
git commit -m "git lfs track"
git push
```

Now commit the bulk of your work:

```bash
git add .
git commit -m "first major commit"
git push
```

You can now visit the `eos9ei3` [GitHub repository](https://github.com/ersilia-os/eos9ei3) and check that your work is publicly available.

### Fetch and serve the model with Ersilia

We are ready to test the model in the context of the Ersilia CLI. To run the model on Caffeine, simply run:

```bash
ersilia fetch eos9ei3
ersilia serve eos9ei3
ersilia api -i "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
```

The input output should lool like this:

```json
{
    "input": {
        "key": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        "input": "Cn1cnc2n(C)c(=O)n(C)c(=O)c12",
        "text": "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
    },
    "output": {
        "outcome": [
            2.297982
        ]
    }
}
```

{% hint style="danger" %}
Debugging the `fetch` and the `api` commands of Ersilia can be very complicated. Please reach out to **@Miquel** directly if you find problems at this step.
{% endhint %}

## TL;DR

In summary, the steps to incorporate a model to the Ersilia Model Hub are the following. We are assuming that the model can be installed in a Conda environment.

1. Download the model from a third party repository to your local machine.
2. Install the model in a dedicated Conda environment and make sure you can run it.
3. Create a new GitHub repository from the `eos-template`. Name this repository with the EOS ID available from the AirTable.
4. Clone the new repository.
5. Place model code in `model/framework` and model parameters in `model/checkpoints`.
6. Write the necessary code to obtain a `run_predict.sh` that simply takes one input file and produces one output file. Be sure to use absolute paths throughout.
7. Edit the `service.py` file, if necessary.
8. Edit the `Dockerfile` file to reflect the installation steps followed in 2.
9. Write the `README` file.
10. Make sure that `.gitattributes` tracks your model parameters.
11. Add, commit and push changes to the model repository.
12. Activate the Ersilia CLI and fetch the model.
13. Serve the model and run the default API.



