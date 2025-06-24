---
description: >-
  This pages provides a deep dive into the structure of the model template for
  new model incorporation.
---

# Model template

## Anatomy of the Ersilia Model Template

Each model in the Ersilia Model Hub is contained within an individual GitHub repository. Each model repository is created using the [**Ersilia Model Template**](https://github.com/ersilia-os/eos-template) upon approval of the Model Request issue. When the new repository is created, please fork it and work on modifying the template from your own user. Open a pull request when the model is ready.

{% hint style="danger" %}
When you have finished the model incorporation, please delete the fork from your own GitHub user. This will prevent abuses of the Git-LFS quota and outdated versions of the models.
{% endhint %}

Below, we describe the main files you will find in the newly created model repository. Note that some of them are automatically updated and you do not have to modify them, like the `README.md`.

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
`eos` identifiers are automatically assigned at repository creation. Please do not modify them.
{% endhint %}

### The `metadata.yml` file

The `metadata.yml` file is where all the model information can be found. This is the only place where you should modify or update the model description, interpretation etc. The Airtable backend, the browsable [interface](https://ersilia.io/model-hub) of the Ersilia Model Hub and the `README.md` file will automatically be updated from the `metadata.yml` upon merge of the Pull Request.&#x20;

YAML fields are constrained by certain parameters. The Pull request triggers a GitHub action that checks the quality of the submitted metadata. If it fails, an explanatory message will be shown on the Action run. Please double check and make amendments if necessary. All accepted and predefined fields in Ersilia Metadata are listed in these .txt [files](https://github.com/ersilia-os/ersilia/tree/master/ersilia/hub/content/metadata). If you wish to include a new one, please open a PR in Ersilia.

**Identifier:** the `eos` identifier described above. It will be automatically filled in. _<mark style="color:red;">Do not modify</mark>._

**Slug:** a one-word or multi-word (linked by a hypen) human-readable identifier, stored as a string, to be used as an alternative to the EOS ID. It will be filled in from the Model Request issue. It can be modified afterwards if necessary.

**Title:** a self-descriptive model title. This field is a single string of minimum 70 characters.

**Description**: minimum information about model type, results and the training dataset. This field is a single string of minimum 200 characters.

{% hint style="info" %}
Some contributors may find it difficult to come up with a good description for the model. You can find some inspiration in [Semantic Scholar](https://semanticscholar.org). This portal provides an AI-based TL;DR short description of many indexed papers.&#x20;
{% endhint %}

**Source:** whether the model runs locally (Local; in your computer if fetched locally, or on an Ersilia cloud for online deployed models) or Online (posts predictions to a server **external** to Ersilia). Note we cannot guarantee privacy of your predictions or input SMILES for Online models. This field is a string with only one accepted value.

**Source type:** where does the model come from. Most typically they are External models (third party) or Internal (developed by the Ersilia team). In some cases we might encounter Replicated models (re-trained using original author's guides). This field is a string with only one accepted value.

**Task**: the ML task performed by the model. This field is typically a list, with one ore more than one entries. The only accepted [tasks](https://github.com/ersilia-os/ersilia/blob/master/ersilia/hub/content/metadata/task.txt) are: `Annotation`, `Representation`, `Sampling`. This field is a string with only one accepted value.

**Subtask:** more granular task description. Annotation models can be `Property calculation or prediction` or `Activity prediction`. Representation models can be `Featurization` or `Projection` and Sampling models can be `Similarity` search or `Generation`. This field is a string with only one accepted value.

**Input**: data type required by the model. Currently all Ersilia models accept only Compound as input. This field is a string with only one accepted value.

**Input dimension:** data format required by the model. Currently all Ersilia models accept only One (in numeric, 1) as input dimension. This does not mean we cannot pass a list of inputs to the model, but the predictions are done on one input, not a combination of them.  This field is an integer with only one accepted value.

**Output:** data type outputed by the model.  The only accepted output formats are: `Boolean`, `Compound`, `Descriptor`, `Distance`, `Experimental value`, `Image`, `Other value`, `Probability`, `Protein`, `Score`, `Text`. This field is a string with only one accepted value.

**Output Dimension:** similar to the input dimension, what is the lenght of the output per each input? . This field is an integer with only one accepted value.

**Output Consistency:** the model produces always the same prediction given the same molecule (Fixed) or not (Variable). Most QSAR models are Fixed, with the exception of Generative models. This field is a string with only one accepted value.

**Interpretation:** provide a brief description of how to interpret the model results. For example, in the case of a binary classification model for antimalarial activity based on experimental IC50, indicate the experimental settings (time of incubation, strain of parasite...) and the selected cut-off for the classification. This field is a single string.

**Tag:** labels to facilitate model search. For example, a model that predicts activity against malaria could have _P.falciparum_ as tag. This field is a list with one or more accepted values since models can have more than one [tag](https://github.com/ersilia-os/ersilia/blob/master/ersilia/hub/content/metadata/tag.txt). Select between one and five relevant from the following categories. tags are predefined in a list at Model Request time. Do not deviate from that. This field is a list of strings.

* Disease: `AIDS`, `Alzheimer`, `Cancer`, `Cardiotoxicity`, `Cytotoxicity`, `COVID19`, `Dengue`, `Malaria`, `Neglected tropical disease`, `Schistosomiasis`, `Tuberculosis`.
* Organism: `A.baumannii`, `E.coli`, `E.faecium`, `HBV`, `HIV`, `Human`, `K.pneumoniae`, `Mouse`, `M.tuberculosis`, `P.aeruginosa`, `P.falciparum`, `Rat`, `Sars-CoV-2`,  `S.aureus`, `ESKAPE`.
* Target: `BACE`, `CYP450`, `GPCR`, `hERG`.&#x20;
* Experiment: `Fraction bound`, `IC50`, `Half-life`, `LogD`, `LogP`, `LogS`, `MIC90`, `Molecular weight`, `Papp`, `pKa`.
* Application: `ADME`, `Antimicrobial activity`, `Antiviral activity`, `Bioactivity profile`, `Lipophilicity`, `Metabolism`, `Microsomal stability`, `Natural product`, `Price`, `Quantum properties`, `Side effects`, `Solubility`, `Synthetic accessibility`, `Target identification`, `Therapeutic indication`, `Toxicity`.
* Dataset: `ChEMBL`, `DrugBank`, `MoleculeNet`, `Tox21`, `ToxCast`, `ZINC`, `TDCommons`.
* Chemoinformatics: `Chemical graph model`, `Chemical language model`, `Chemical notation`, `Chemical synthesis`, `Compound generation`, `Descriptor`, `Drug-likeness`, `Embedding`, `Fingerprint`, `Similarity`.&#x20;

**Biomedical area:** the pertinent area of research or disease targeted. More than one can be selected from the same model, for example ADMET and Malaria, if both are relevant. Any can be used when the model could be applied to all fields of disease (for example a featuriser). This field is a list of strings.

**Target organism:** if existing, the pathogen the model is related to, for example _Plasmodium falciparum,_ or organism in case of ADMET, for example _Mus musculus_ or _Homo sapiens._ If the model is unrelated to any organism, select Not applicable. This field is a list of strings.

**Publication type:** peer-reviewed, preprint or other. Most Ersilia models are peer-reviewed, and only on exception "other" is accepted. This field is a string with only one accepted value.

**Publication year:** year of publication of the original model. This field is an integer with only one accepted value.

**Publication:** link to the original publication. Please refer to the journal page whenever possible, instead of PubMed, ResearchGate or other secondary webs. If the model is not published in a scientific journal, please contact Ersilia maintaners. This field is a string with URL format.

**Source Code:** link to the original code repository of the model. If the model is a de-novo incorporation without another source code, please link the ML package used to train the model. This field is a string with URL format.

**License:** License of the original code. Only SPDX abreviations are accepted. If the code was released without a license, please add `None` in this field. Make sure to abide by requirements of the original license when re-licensing or sub-licensing third-party author code (such as adding the license file together with the original code). This field is a string with only one accepted value.

{% hint style="info" %}
If the predetermined fields are not sufficient for your use case, you can open a pull request to include new ones to our [repository](https://github.com/ersilia-os/ersilia/tree/master/ersilia/hub/content/metadata). Please do so only if strictly necessary (for example, if a disease is not already in the Tag field).

Ersilia maintainers will review and approve / reject PRs for additions to the existing lists of approved items.
{% endhint %}

{% hint style="danger" %}
Note that these fields are filled in as Python strings, therefore misspellings or lower / uppercases will affect their recognition as valid values.
{% endhint %}

### The [`README`](https://github.com/ersilia-os/eos-template/blob/main/README.md) file

The `README.md` file is where we give basic information about the model. It reads from the metadata.yml file and it will be automatically updated thanks to a GitHub Action once the Pull Request is approved.

Please do not modify it manually.

### The [`LICENSE`](https://github.com/ersilia-os/eos-template/blob/main/LICENSE) file

By default, all code written in contribution to Ersilia should be licensed under a [GPLv3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). The main `LICENSE` file of the repository, therefore, will be a GPLv3 as specified by [GitHub.](https://docs.github.com/en/communities/setting-up-your-project-for-healthy-contributions/adding-a-license-to-a-repository)

However, the license notices for code developed by **third parties** must be kept in the respective folders where the third-party code is found.

### The [`install.yml`](https://github.com/ersilia-os/eos-template/blob/main/install.yml) file

Ersilia uses an `install.yml` file to specify installation instructions. The YAML syntax is easy to read and maintain. This file specifies which Python version to use to build a conda environment, or a Docker image for the model.&#x20;

This dependency configuration file has two top level keys, namely, `python`, and `commands.` They dependencies are to be specified in the following manner:

* `python` key expects a string value denoting a python version (eg `"3.10"`)
* `commands` key expects a list of values, each of which is a list on its own, denoting the dependencies required by the model. Currently, dependencies  from`pip` and `conda` are supported.
* `pip` dependencies are expected to be three element lists in the format `["pip", "library", "version"]`
* `conda` dependencies are expected to be four element lists in the format `["conda", "library", "version", "channel"]`, where channel is the conda channel to install the required library.
* When the model is run from source, Ersilia always defaults to creating a conda environment for the model to provide isolation. However, when the model is Dockerized, whether conda is used in that process depends entirely on there being conda dependencies in this file.&#x20;

The `install.yml` available in the Ersilia Model Template only defines two top-level keys as follows:

```yaml
python: "3.10"        # REQUIRED; value must be a string
commands:            # REQUIRED; value must be a sequence (list)
  - ["pip", "package", "1.2.3"]
  - ["pip", "git+https://…"]
  - ["conda", "install", "-c", "channel", "pkg=4.5.6"]
```

The commands key only allows the following entries:

* String: any shell command passed verbatim

```yaml
- "apt-get update"
- "mkdir -p /data"
```

* List with pip: a two-element list for Git URLs that begin with `git+https://`, `git+ssh://`, or `https://` or a list with three or more elements containing the package name (compulsory), package version  (compulsory) and optional pip flags

```yaml
- ["pip", "git+https://github.com/user/repo.git"]
- ["pip", "mypkg", "0.1.2", "--upgrade", "--no-deps"]
```

* List with conda: users can specify a long version (recommended) where channels are indicated by `-c` or use a shorthand form with the following order: `conda`, package name, version string, then any channels.  If the `-y` flag is not added, it will be automatically appended in both cases. Packages must have the version specified either as `package=version` or as `"[package", "version"].` If no channel is specified, default will be used

```yaml
- ["conda", "install", "-y", "-c", "channel", "package=version"] #long one channel
- ["conda", "install", "-y", "-c", "channel1", "-c", "channel2","package=version"] #long two channels
- ["conda", "package", "version", "channel1", "channel2", "-y"] #short version with channels specified
- ["conda", "package", "2024.3.6"] # short version no channel specified
```

This example demonstrates an installation instructions for an environment using Python 3.10.

{% code title="install.yml" %}
```yaml
python: "3.10"
commands:
    - ["pip", "rdkit-pypi", "2022.3.1b1"]
    - ["conda", "pandas", "1.3.5", "default"]
```
{% endcode %}

In this example, the `rdkit-pypi==2022.3.1b1` will be installed using `pip`, while `pandas==1.3.5` will be installed using `conda` through the default package channel on conda.

The `install.yml` file can contain as many commands as necessary. Please limit the packages to the bare minimum required, sometimes models have additional packages for extra functionalities that are not required to actually run the model. It is good practice to trim to the minimum the package dependencies to avoid conflicts. Always pin the version of the package for model persistency.

{% hint style="warning" %}
The `install.yml` file contains the installation instructions of the model. Therefore, the content of this file can be very variable, since each model will have its own dependencies.
{% endhint %}

### The [`model`](https://github.com/ersilia-os/eos-template/tree/main/model) folder

The `model` folder is the most important one. It contains two sub-folders:

* `framework`: In this folder, we keep all the necessary code to run the model (assuming dependencies are already installed).
* `checkpoints`: In this folder, we store the model data (pretrained model parameters, scaling data, etc).

{% hint style="danger" %}
The `model` folder **should not** contain anything other than the `framework` and `checkpoints` subfolder. When the Ersilia CLI eventually fetches the model, it does a reorganization of the code and the only subfolders it keeps are these two. Any other file or folder at the `model/` directory level will be overlooked.
{% endhint %}

{% hint style="info" %}
Often, the separation between `framework` and `checkpoints` is not easy to determine. Sometimes, models obtained from third parties have model data embedded within the code or as part of the repository. If those are impossible to separate, keep the data in the `framework` subfolder, and leave the `checkpoints` subfolder empty.
{% endhint %}

The `framework` subfolder contains at least one Bash file, named `run.sh`. This file will run as follows:

```bash
bash run.sh [FRAMEWORK_DIR] [DATA_FILE] [OUTPUT_FILE]
```

Unless strictly necessary, the `run.sh` file should accept three and only three arguments, namely `FRAMEWORK_DIR`, `DATA_FILE` and `OUTPUT_FILE`. In the current template, we provide the following example:

{% code title="run.sh" %}
```bash
python $1/code/main.py -i $2 -o $3
```
{% endcode %}

In this case, a Python file located in the `[FRAMEWORK_DIR]/code` folder is executed, taking as input (`-i`) the `DATA_FILE` and giving as output (`-o`) the `OUTPUT_FILE`.

To understand this further, we now need to inspect the step `main.py` file in the step above, in more detail. The current template proposes the following script:&#x20;

{% code title="code/main.py" %}
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

Most of the work of the model contributor will be to work on this or similar scripts. In the template, we provide a dummy model (i.e. add a fixed value to the molecular weight). This dummy model can can be already defined within the script (`my_model`). However, in real world cases, the model will most likely be loaded from a third party Python library, or from a script placed in the same directory.

To summarize, in the template, we provide a structure that follows this logic:

1. The `run.sh` script executes the Python `main.py` script.
2. The `main.py` script:
   * Defines the model code.
   * Loads parameters from `checkpoints`.
   * Reads an input file containing SMILES (with header).
   * Runs a model that calculates molecular weight and adds an integer defined by the parameters.
   * Writes an output file containing one column corresponding to the output value (with a header).

{% hint style="info" %}
In the template, the example provided is very simple. Depending on the model being incorporated, the logic may be different. For example, many third party models already contain a command-line option, with a specific syntax. In these cases, you may want to write scripts to adapt the input and the output, and then execute the model as-is.

Each script will be one `.py` file, we can create as many as necessary and rename them appropriately (see below for examples).
{% endhint %}

### The `.gitattributes` file

We use Git LFS to store large files (over 100 MB). Typically, these files are model parameters. Files to be stored in Git LFS should be specified in the `.gitattributes` file. Examples for Git LFS tracking of `csv`, `h5`, `joblib`, `pkl`, `pt` and `tsv` are:

```
*.csv filter=lfs diff=lfs merge=lfs -text
*.h5 filter=lfs diff=lfs merge=lfs -text
*.joblib filter=lfs diff=lfs merge=lfs -text
*.pkl filter=lfs diff=lfs merge=lfs -text
```

{% hint style="warning" %}
To allow model contributors to push Git LFS files to Ersilia's repository (counts towards Ersilia's quota) a mock.csv file will be creatd and tracked. Make sure to eliminate the mock.csv file and, if not needed for other files, the .gitattributes as well.
{% endhint %}

## The Legacy Model Template

At the time of writing this tutorial, the Ersilia Model Hub has approximately 150 models developed with this template. As mentioned above, this legacy template served as an inspiration for the current workflow, however there are several key differences that should be called out. These differences largely pertain to how dependencies are specified in these models, the tools used to create a model server, and source files facilitating that process.

### The [`Dockerfile`](https://github.com/ersilia-os/eos-template/blob/42ce4063e67122968c3c948bd8ea142ac621c105/Dockerfile) file

Ersilia uses a `Dockerfile` file to specify installation instructions. The reason for this is that Docker provides the maximum level of isolation possible (i.e. a container), which may be needed to run models in some systems. However, in most practical scenarios, a Docker container will not be necessary and a Conda environment, or even a Virtualenv environment, will suffice. The Ersilia CLI will decide which isolation level to provide depending on the content of the `Dockerfile:`

The `Dockerfile` uses a templated file, where the lines beginning with `RUN` are scanned; installer commands are normalized, all other shell commands are preserved.

1. Non-`pip`/`conda` RUN commands\
   – Entire line (after `RUN`) is retained as a single string.
2. Pip installs:
   * Token count must be ≥3 (`pip`, `install`, `pkg_spec`)
   * `pkg_spec` starting with `git+` → `["pip", "git+…"]`
   *   Otherwise, `pkg_spec` must contain `==version`; its absence will raise the following error:

       ```
       ValueError("pip install must specify version or git URL")
       ```
   * Additional tokens → appended as flags, resulting in\
     `["pip", "pkg", "ver", …flags]`
3. Conda installs:

* Tokens must begin with `conda install` and package version must be specified
* Entire token list returned, e.g.:\
  `["conda", "install", "-y", "numpy=1.23.0", "-c", "conda-forge"]`

The python version is specified from the first line, even though BentoML is no longer in use.

<pre class="language-docker"><code class="lang-docker"><strong>Dockerfile 
</strong>FROM bentoml/model-server:0.11.0-py310
…

RUN pip install mypkg==0.3.1
RUN pip install mypkg[option]==0.3.1
RUN pip install git+https://github.com/user/repo.git
RUN conda install -y numpy=1.23.0 -c conda-forge -c other-channel -y

WORKDIR /repo 
COPY . /repo
</code></pre>

### The [`service`](https://github.com/ersilia-os/eos-template/blob/42ce4063e67122968c3c948bd8ea142ac621c105/src/service.py) file

The service file is located in `src/service.py`. It contains the necessary code to facilitate model bundling with BentoML.

There are three main classes in the `service` file, namely `Model`, `Artifact` and `Service`.

#### The `Model` class

This class is simply a wrapper for the AI/ML model. Typically, when incorporating **external** (type 1) models, the `run.sh` script will already capture the logic within the `Model` class, in which case the `Model` class is simply redundant. However, when incorporating **internally developed** (types 2 and 3) models into the hub, we can make use of the artifacts for standard modeling frameworks (e.g. sklearn, PyTorch, and Keras) provided by BentoML, and the `Model` class becomes necessary for BentoML compatibility. Hence, the `Model` class enables generalization between these types of model additions.

Typically, the central method of the `Model` class is the `run` method.

```python
class Model (object):
    ...
    def run(self, smiles_list):
        ...
```

In this case, the model takes as input a list of molecules represented as SMILES strings. This is the standard input type for models focused on chemistry data as input.

{% hint style="warning" %}
Models incorporated with this workflow do not allow for multiple endpoints, or "methods". All models developed in this manner only have a single `run` API.
{% endhint %}

In its simplest form, the `Model` class just points Ersilia to the `model` directory and then creates a Bash file to execute the necessary commands to run the model. It is actually a very simple class, although it may look overwhelming at first. We break it down below:

First, a temporary directory is created:

```python
 class Model(object):
    ...
    def run(self, smiles_list):
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
    def run(self, smiles_list):
        ...
        with open(data_file, "w") as f:
            f.write("smiles"+os.linesep)
            for smiles in smiles_list:
                f.write(smiles+os.linesep)
        ...
```

Now we already have the input file of the `run.sh`script, located in the `model/framework/` directory, as specified [above](model-template.md#the-model-folder). The following creates a dummy Bash script in the temporary directory and runs the command from there. The output is saved in the temporary directory too. [Remember](model-template.md#the-model-folder) that the `run.sh` script expects three arguments, `FRAMEWORK_DIR`_,_ `DATA_FILE` and `OUTPUT_FILE`.

```python
class Model(object):
    ...
    def run(self, smiles_list):
        ...
        run_file = os.path.join(tmp_folder, self.RUN_FILE)
        with open(run_file, "w") as f:
            lines = [
                "bash {0}/run.sh {0} {1} {2}".format(
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
    def run(self, smiles_list):
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

You will see that, in the template, pointers to potential edits are highlighted with the tag `# EDIT` . Necessary edits relate to the the format of the input data, or the serialization to JSON format from the output data.

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

You <mark style="color:red;">**don't have to modify**</mark> this class.

#### The `Service` class

This class is used to create the service. The service exposes the `run` API:&#x20;

```python
@artifacts([Artifact("model")])
class Service(BentoService):
    @api(input=JsonInput(), batch=True)
    def run(self, input: List[JsonSerializable]):
        input = input[0]
        smiles_list = [inp["input"] for inp in input]
        output = self.artifacts.model.run(smiles_list)
        return [output]
```

By default, Ersilia works with JSON inputs, which are deserialized as a SMILES list inside the API, in this case. The deafult API is `run`. The general rule is, <mark style="color:red;">**do not modify**</mark> it.
