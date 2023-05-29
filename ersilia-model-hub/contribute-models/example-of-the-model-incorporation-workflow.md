# Model incorporation workflow

## Overview

The model incorporation workflow is streamlined by a series of GitHub Actions. Each push to Ersilia's codebase triggers a series of events that automatically maintain our platform updated across the different services we use to maintain our infrastructure (GitHub, Airtable, DockerHub and AWS).

There are two sets of workflows. The ones triggered by changes in the Ersilia Model Hub itself, and the ones triggered by each individual model. Here, we will focus on the workflow triggered when a new model is being incorporated:&#x20;

<figure><img src="../../.gitbook/assets/workflow (1).png" alt=""><figcaption><p>Diagram of the Ersilia Workflow at new model incorporation. Purple indicates automated GitHub Actions, bold indicates Ersilia team and contributor actions.</p></figcaption></figure>

## A toy example

If you want to test the workflow with a toy example, we have prepared an [Ersilia Demo Repository](https://github.com/ersilia-os/eos-demo). This repository is focused on illustrating the GitHub Actions workflows to be followed in order to successfully incorporate a model to the Ersilia Model Hub.

### 1. Create a model request at Ersilia

1. Go to the Ersilia main repository [issues page](https://github.com/ersilia-os/ersilia/issues).
2. Click on **New issue**. Then 🦠 **Model Request** (**Get started**)
3. For the purpose of this demo, you can use the following information:
   * 🦠 **Model Request**: Demo Malaria Model
   * **Model Name**: Demo Malaria Model
   * **Model Description**: Prediction of the antimalarial potential of small molecules. This model was originally trained on proprietary data from various sources, up to a total of >7M compounds. The training sets belong to Evotec, Johns Hopkins, MRCT, MMV - St. Jude, AZ, GSK, and St. Jude Vendor Library. In this implementation, we have used a teacher-student approach to train a surrogate model based on ChEMBL data (2M molecules) to provide a lite downloadable version of the original MAIP.
   * **Slug**: demo-malaria-model
   * **Tag**: Malaria,P.falciparum
   * **Publication**: [https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00487-2](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00487-2)
   * **Code**: [https://www.ebi.ac.uk/chembl/maip/](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00487-2)
   * **License**: GPL-3.0

### 2. Wait until model approval

1. The Ersilia team will revise your model requests and likely start a public discussion around it.
2. At some point, your model will be approved. You will see an `/approve` comment in the GitHub issues.
3. Approval will trigger some GitHub Actions. Eventually, the `ersilia-bot` will post an informative message in your issue. Importantly, this message will contain a link to a new model repository placeholder, named, for example, `ersilia-os/eosXabc`. This code is arbitrarily assigned by Ersilia and will change every time. You should not worry about creating one manually and you should never modify it.

### 3. Fork the model repository

1. Go to the model repository page: https://github.com/ersilia-os/eosXabc, in this case.
2. Fork the repository to your username.
3. Clone the forked repository. This should create a `eosXabc` folder in your local filesystem.

### 4. Use this demo model

1. For this demo, you have to clone the current repository (`ersilia-os/eos-demo`). This will create an `eos-demo` folder in your local filesystem.
2. Run the following script to populate your forked model with the demo data: `python /path/to/eos-demo/populate.py /path/to/eosXabc`.
3. Your `eosXabc` has been populated with model parameters, dependencies and extra metadata. It is now ready for commit.

### 5. Make a pull request to Ersilia

1. Commit changes and push changes to `eosXabc`.
2. Open a PR to the `main` branch at `ersilia-os/eosXabc`. GitHub Actions workflows will be triggered to ensure that your code works as expected.

### 6. Wait until model is merged

1. The Ersilia team will revise your PR and merge it eventually. More GitHub Actions workflows will be triggered at this point.
2. Once the model is merged, you should see it in [Ersilia's AirTable](https://airtable.com/shrNc3sTtTA3QeEZu).

### 7. Assist with curation and publication

1. The `ersilia-bot` will open a new issue at `ersilia-os/eosXabc`. As you will see, someone from the Ersilia community will be assigned as a reviewer of the model.
2. If you are a member of the [Ersilia Slack workspace](https://ersilia-workspace.slack.com/), then you may also see activity triggered around your model.

## A real-world example

Now that we have an idea of the contents of the [Ersilia Model Template](https://github.com/ersilia-os/eos-template), we will follow the example of a simple but widely used model to calculate the **synthetic accessibility** of small molecule compounds. Synthetic accessibility measures the feasibility of synthesizing a molecule in the laboratory. In [2009, Peter Ertl presented the synthetic accessiblity (SA) score](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8), based on measures of molecular complexity and occurrence of certain fragments in the small molecule structure. High (greater than 6) SA scores denote difficult-to-synthesize molecules, and low (lower than 3) SA scores suggest that the molecule will be easy to synthesize.

### 1. Open a Model Request Issue

Please fill in the[ issue](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&template=model\_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) fields as accurately as possible and wait for review and approval by one of the Ersilia maintainers.

#### Read the publication

It is important that you read the original publication in order to understand the training data, the limitations of the model, the validation performed and the characteristic of the algorithm, among other details. In [this case](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8), this is a classic (old) publication from a Novartis team. They analyzed small fragments in the PubChem database and devised a molecular complexity score that takes into account molecule size, presence of non-standard structural features, presence of large rings, etc. The authors validated the model by comparing their SA score with synthetic feasibility as estimated by medicinal chemists for a set of 40 molecules.

#### Find model code and parameters

Code to calculate the SA score does not seem to be available from the publication. Fortunately, though, the RDKit library, in its contributions module, contains an implementation of the SA score. The code can be found [here](https://github.com/rdkit/rdkit/tree/master/Contrib/SA\_Score). This RDKit-based implementation was developed in 2013 by Peter Ertl and Greg Laundrum.

{% hint style="success" %}
Both the link to the code and to the original publication are accessible from the Ersilia Model Hub AirTable database.
{% endhint %}

### 2. Run the code outside Ersilia

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

### 3. Fork the new model repository based on the Ersilia Model Template

Now that we know that the code can run in our local machine, we can fork the new repository created by the <mark style="color:green;">Model Request Workflow</mark> and start working on it.

#### Fork and Clone the new repository

Always fork the repository to your user, and clone it from there in our local machine. Let's do it in the `~/Desktop`:

```bash
cd ~/Desktop
git clone https://github.com/user-github/eos9ei3.git
cd eos9ei3
```

#### Migrate code and parameters

Let's now place the code and the parameters in the `model` folder (in the `framework` and `checkpoints` sub folders, respectively):

```bash
cd ~/Desktop
cp SA_Score/sascorer.py eos9ei3/model/framework/code/.
cp SA_Score/fpscores.pkl.gz eos9ei3/model/checkpoints/.
```

{% hint style="info" %}
The checkpoints contains a dummy `checkpoints.joblib` that can be deleted.
{% endhint %}

{% hint style="danger" %}
Note that here we are migrating code and parameters to different folders. This may cause critical errors if code expects to find parameters at a certain relative location. Try to locate the pointers to the model parameters and change the paths. Only, and exceptionally, if the model architecture is too complex we can keep code and parameters in the `/framework` folder. Please ask for permission to Ersilia's team before doing it.
{% endhint %}

#### Write framework code

Now it is time to write some code. Here we will follow the description of the `model` folder given in the [Model Template page](model-template.md#the-model-folder).

The `eos-template` provides a `main.py` that will guide us through the deployment steps. We will adapt the `main.py` to our specific needs, and run it from the `run.sh` file present in the /model directory. The arguments are, respectively, the path to the file, the input and the output. Ersilia takes care of passing the right arguments to the `run.sh` file. Please do not modify it.

```bash
python $1/code/main.py $2 $3
```

Going back to our model of interest, we have identified three necessary steps to run the model:

1. Input adapter
2. SAScorer (in this case, the `sascorer.py` which is already written)
3. Output adapter

**Write the input adapter**

By default, for chemical compound inputs, Ersilia uses single-column files with a header (see the `service.py` [file](model-template.md#the-service-file)). However, the `sascorer.py` [expects](example-of-the-model-incorporation-workflow.md#test-the-model) a two-column file. Let's write an **input** **adapter**:

{% code title="input_adapter.py" %}
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

The script creates an intermediate `tmp_input.smi` file that can be used as input for `sascorer.py`. We can keep this as a separate file under /code, but since it is a small function, we will write it inside the `main.py` file:

{% code title="code/main.py from line 20" %}
```python
# read SMILES from .csv file, assuming one column with header
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

**Make sure that parameters are read**

So far, we haven't pointed to the model parameters. When [migrating code and parameters](example-of-the-model-incorporation-workflow.md#migrate-code-and-parameters), we separated the `sascore.py` file and the `fpscores.pkl.gz` file.

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

**Run the model inside main.py**

We now have the input adapter and the model code and parameters. We simply need to run the model in main.py by calling the `sascorer.py.`

We first make sure to import the necessary packages and delete the non-necessayr (in this case, the MW by rdkit). The sascorer uses time to measue how long did the model take, as well as the functions readFragmentScores and processMols.

The input file is no longer passed as a sys.argv, so we modify it with the temporal input file we just created

{% code title="code/main.py" %}
```python
# imports
import os
import csv
import sys
import time
from rdkit import Chem
from sascorer import readFragmentScores, processMols

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# read SMILES from .csv file, assuming one column with header
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

# run model
t1 = time.time()
readFragmentScores("fpscores")
t2 = time.time()

suppl = Chem.SmilesMolSupplier("tmp_input.smi")
t3 = time.time()
R = processMols(suppl)
t4 = time.time()

print('Reading took %.2f seconds. Calculating took %.2f seconds' % ((t2 - t1), (t4 - t3)),
        file=sys.stderr)


```
{% endcode %}

The `processMols()` function was simply printing the output without saving the results, so we have modified it in the `sascorer.py` to get the output:

```python
def processMols(mols):
    print('smiles\tName\tsa_score')
    R = []
    for i, m in enumerate(mols):
        if m is None:
            continue
        s = calculateScore(m)
        smiles = Chem.MolToSmiles(m)
        print(smiles + "\t" + m.GetProp('_Name') + "\t%3f" % s)
        R += [[smiles, m.GetProp("_Name"), "\t%3f" % s]]
    return R
```

#### Write the output adapters

We need to understand the output of the model in order to collect it correctly. The easiest will be to add a print(R) statement in main.py and run the `run.sh` file with a mock test file.

```bash
cd model/framework
python run.sh . ~/test.csv ~/output.csv
```

We observe that the output provided by `sascorer.py` has three columns (tab-separated), so we need to adapt it.

<table><thead><tr><th width="462.3333333333333">smiles</th><th>Name</th><th>sa_score</th></tr></thead><tbody><tr><td>C(F)Oc1ccc(-c2nnc3cncc(Oc4ccc5ccsc5c4)n23)cc1</td><td>mol0</td><td>2.823995</td></tr><tr><td>C(F)Oc1ccc(-c2nnc3cncc(OCC[C]4BBBBBBBBBB[CH]4)n23)cc1</td><td>mol1</td><td>5.757383</td></tr><tr><td>Cn1cc2ccc(Oc3cncc4nnc(-c5ccc(OC(F)F)cc5)n34)cc2n1</td><td>mol2</td><td>2.910502</td></tr></tbody></table>

Likewise, we can subsitute this piece of code in `main.py`:

{% code title="code/main.py  template line29" %}
```python
outputs = []
for r in R:
    outputs += [r[-1].strip()]
        
#check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["sa_score"])  # header
    for o in outputs:
        writer.writerow([o])
```
{% endcode %}

Finally, let's add one more line to main.py to clean up the temporary files we have created:

```python
os.remove("tmp_input.smi")
```

### 4. Run the model locally

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

To test the model, we simply have to execute the `run.sh` script.

```
cd ~/Desktop
FRAMEWORK_PATH="eos9ei3/model/framework/"
bash $FRAMEWORK_PATH/run.sh $FRAMEWORK_PATH molecules.csv output.csv
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

#### Edit the `service.py` file, if necessary

The `service.py` file provided by default in the template manages chemistry inputs and expects tabular outputs. Therefore, in principle, you do not have to modify this file.

{% hint style="info" %}
Modifying the `service.py` file is intended for advanced users only. Please use the Slack `#internships` channel if you think your model of interest requires modification of this file.
{% endhint %}

#### Edit the `Dockerfile` file

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

#### Write the `metadata.json` file

Don't forget to document the model. Read the [instructions to write the `metadata` file](example-of-the-model-incorporation-workflow.md#the-metadata.json-file) page. Feel free to ask for help in the Slack `#internships` channel.

The metadata.json for this model should read like:&#x20;

```json
{
    "Identifier": "eos9ei3",
    "Slug": "sa-score",
    "Status": "Ready",
    "Title": "Synthetic accessibility score",
    "Description": "Estimation of synthetic accessibility score (SAScore) of drug-like molecules based on molecular complexity and fragment contributions. The fragment contributions are based on a 1M sample from PubChem and the molecular complexity is based on the presence/absence of non-standard structural features. It has been validated comparing the SAScore and the estimates of medicinal chemist experts for 40 molecules (r2 = 0.89). The SAScore has been contributed to the RDKit Package.\n",
    "Mode": "Pretrained",
    "Input": [
        "Compound"
    ],
    "Input Shape": "Single",
    "Task": [
        "Regression"
    ],
    "Output": [
        "Score"
    ],
    "Output Type": [
        "Float"
    ],
    "Output Shape": "Single",
    "Interpretation": "Low scores indicate higher synthetic accessibility",
    "Tag": [
        "Synthetic accessibility",
        "Chemical synthesis"
    ],
    "Publication": "https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8",
    "Source Code": "https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score",
    "License": "BSD-3.0",
    "Contributor": "github-user"
```

### 5. Run the local model inside Ersilia

Before committing our new model to Ersilia, we must check it will work within the Ersilia environment. To do so, we have a very convenient option at model fetch time, `--repo_path` that allows us to specify a local path to the model we are fetching (so, instead of looking for it online it will use the local folder). <mark style="color:purple;">**It is crucial to complete this step**</mark> before committing the model to GitHub.

```
conda activate ersilia
ersilia -v fetch eos9ei3 --repo_path ~/Desktop/eos9ei3
ersilia serve eos9ei3
ersilia predict -i molecules.csv -o output.csv
```

If this runs without issues, the model is ready to be incorporated. If not, please go back to step 4 and revise the model indeed is running without issues.

### 6. Commit changes to the repository

We are now ready to commit changes, first to our fork, and then to the main Ersilia repository by opening a pull request. Before doing so, complete the steps below:

#### Cleanup mock files

Probably, there is a few files, such as `mock.csv`, that are no longer needed (this is used solely to initialize Git-LFS). Please remove them before committing and also eliminate the Git-LFS tracking from `.gitattributes` (see below for more on Git-LFS)

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

#### Open a pull request

Once the model is ready, open a pull request to merge your changes back into the main repository. This will trigger a series of GitHub Actions:

* <mark style="color:green;">Security workflow</mark>: makes sure that no private key is released with the model
* <mark style="color:green;">Json syntax check:</mark> controls that the metadata file does not contain Json syntax errors (does not check the content, only syntax)
* <mark style="color:green;">Model Test on PR</mark>: this workflow will first test that the metadata.json has the correct fields (if it doesn't, it will fail. Please look at the Action Run to get more information on why it has failed). If the metadata is correct, it will then go onto installing Ersilia and testing the model.

If the Actions at Pull request fail, please check them and work on debugging them before making a new pull request. <mark style="color:purple;">**Ersilia maintainers will only merge PR's that have passed all the checks.**</mark> Once the PR has the three green checks, the PR will be merged. This triggers a <mark style="color:green;">Model Test on Push</mark> action that will:

* Test the model once more
* Update the README.md and AirTable metadata&#x20;
* Open a new issue requesting two Ersilia Community members to test the model.

This workflow is also triggered each time there is a push to the repository, to ensure changes to the code do maintain model functionality and metadata is not outdated. If the model testing works, two final actions will be triggered:

* <mark style="color:green;">Upload model to Dockerhub</mark>: the model will be packaged in a docker image and made available via the [Ersilia DockerHub](https://hub.docker.com/orgs/ersiliaos) page
* <mark style="color:green;">Upload model to S3</mark>: the model is zipped and uploaded to S3, to facilitate upload and download from the CLI and avoid incurring Git-LFS bandwith problems.

You can now visit the `eos9ei3` [GitHub repository](https://github.com/ersilia-os/eos9ei3) and check that your work is publicly available.

### 7. Fetch and serve the model with Ersilia

We are ready to test the model in the context of the Ersilia CLI. To run the model, simply run:

```bash
ersilia fetch eos9ei3
ersilia serve eos9ei3
ersilia api -i "Cn1cnc2n(C)c(=O)n(C)c(=O)c12"
```

The input output should look like this:

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

As mentioned, the workflow will also trigger a request for model testing to members of the Ersilia community via a GitHub issue in the same repository. The original model contributor should make sure the modle is working for different users and answer any questions or issues that might arise during model testing. If amends must be made, the original model contributor should work on those

#### 8. Clean up

Please, help us keep a healthy environment and avoid duplication of files and consuming of our Git LFS quota. Delete the repository fork after the model has been successfully incorporated.

## TL;DR

In summary, the steps to incorporate a model to the Ersilia Model Hub are the following. We are assuming that the model can be installed in a Conda environment.

* Download the model from a third party repository to your local machine.
* Install the model in a dedicated Conda environment and make sure you can run it.
* Open a Model Request issue in the Ersilia Model Hub code repository. Wait for approval to automatically obtain a new model repository from the `eos-template`
* Fork the new repository.
* Place model code in `model/framework` and model parameters in `model/checkpoints`.
* Write the necessary code to obtain a `run.sh` that simply takes one input file and produces one output file. Be sure to use absolute paths throughout.
* Edit the `service.py` file, if necessary.
* Edit the `Dockerfile` file to reflect the installation steps followed in 2.
* Update the `metadata.json` following the guidelines
* Make sure that `.gitattributes` tracks your model parameters.
* Open a PR to the the model repository, and check that all tests are passed. If not, try to identify the bug to solve it.
* Once the PR passes all the tests and is merged, activate the Ersilia CLI and fetch the model.
* Serve the model and run the default API.
* Check other users can also run the model.
* Delete your fork.
