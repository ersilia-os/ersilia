---
description: Step by step for the incorporation of models in the Ersilia Model Hub.
---

# Model incorporation workflow

{% hint style="info" %}
This page is **under construction** :construction\_worker:.
{% endhint %}

## A toy example

This example has been extracted from the [Ersilia Demo Repository](https://github.com/ersilia-os/eos-demo). This repository is focused on illustrating the GitHub Actions workflows to be followed in order to successfully incorporate a model to the Ersilia Model Hub.

### Steps

#### Create a model request at Ersilia

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

#### Wait until model approval

1. The Ersilia team will revise your model requests and likely start a public discussion around it.
2. At some point, your model will be approved. You will see an `/approve` comment in the GitHub issues.
3. Approval will trigger some GitHub Actions. Eventually, the `ersilia-bot` will post an informative message in your issue. Importantly, this message will contain a link to a new model repository placeholder, named, for example, `ersilia-os/eosXabc`. This code is arbitrarily assigned by Ersilia and will change every time. You should not worry about creating one manually and you should never modify it.

#### Fork the model repository

1. Go to the model repository page: https://github.com/ersilia-os/eosXabc, in this case.
2. Fork the repository to your username.
3. Clone the forked repository. This should create a `eosXabc` folder in your local filesystem.

#### Use this demo model

1. For this demo, you have to clone the current repository (`ersilia-os/eos-demo`). This will create an `eos-demo` folder in your local filesystem.
2. Run the following script to populate your forked model with the demo data: `python /path/to/eos-demo/populate.py /path/to/eosXabc`.
3. Your `eosXabc` has been populated with model parameters, dependencies and extra metadata. It is now ready for commit.

#### Make a pull request to Ersilia

1. Commit changes and push changes to `eosXabc`.
2. Open a PR to the `main` branch at `ersilia-os/eosXabc`. GitHub Actions workflows will be triggered to ensure that your code works as expected.

#### Wait until model is merged

1. The Ersilia team will revise your PR and merge it eventually. More GitHub Actions workflows will be triggered at this point.
2. Once the model is merged, you should see it in [Ersilia's AirTable](https://airtable.com/shrNc3sTtTA3QeEZu).

#### Assist with curation and publication

1. The `ersilia-bot` will open a new issue at `ersilia-os/eosXabc`. As you will see, someone from the Ersilia community will be assigned as a reviewer of the model.
2. If you are a member of the [Ersilia Slack workspace](https://ersilia-workspace.slack.com/), then you may also see activity triggered around your model.

## A real-world example

## **Model example: Prediction and Optimization of Human Plasma Protein Binding (PPB) of Compounds**

### **1.  Choose the model to incorporate into the Ersilia model Hub.**

Currently, Ersilia has a database of prospective models from the scientific literature developed by third parties, which are hosted on the Airtable.

We can choose a model to work from the [Ersilia airtable model](https://airtable.com/shrTpe45mLKqaHXsc/tblvRnXT57xljIjr0) and can filter by highest priority and choose a model of greater interest.

If, on the other hand, you want to request a new model that is not in the previous database, you can request the model through this link: [https://airtable.com/shroQLlkcmDcC0xzm](https://airtable.com/shroQLlkcmDcC0xzm), or contact the organizers of Ersilia (Gemma and Miquel).

Next, we will use as an example the [IDL-PPBopt ](https://airtable.com/shrTpe45mLKqaHXsc/tblvRnXT57xljIjr0/viwn2g36BdOYHqstG/recguCJqFN1AQzVf0)model, which was chosen from the Ersilia model hub database in Airtable. There we can find the links to the publication and the links to source code of the model.

### **2.  Model Interpretation**

For the interpretation of the model, it is very important to carefully read the original publication of the model since it is from here that we can see model viability for Ersilia. It is important at this point to identify several relevant characteristics such as:

* Understand the general purpose of the model. Reflect on the importance of the model for Ersilia and if this is aligned with its main objective.
* Recognize the type of Input.&#x20;
* Recognize the type of model: Check if the model is pre-trained, retrained, in-house, or Online.
  * In this case check if the author is providing the checkpoints, or if the author is not providing the pre-trained parameters then check if he provides the data set to retrain the model.
* Identify if the model is classification, regression, representation, generative, etc.
* Type of the Output.&#x20;
* Model Code Identification:
  * Identification of the model requirements (packages to install and python versions compatible with Ersilia).

For this [IDL-PPBopt ](https://pubs.acs.org/doi/10.1021/acs.jcim.2c00297)model, the objective is to predict the plasma protein binding ([PPB](https://en.wikipedia.org/wiki/Plasma\_protein\_binding) \*) property based on an interpretable deep learning method.

When the drug enters the body, it interacts with plasma proteins, being able to remain bound to the plasma protein or remain free and find its pharmacological target. The plasma protein-bound form is released slowly, prolonging the duration of action of the drug.

In this model to predict the Plasma Protein Binding property. The author first creates a deep learning model interpretable using the attentive fingerprinting algorithm (AFP) to predict the fraction of PPB (this fraction measures the binding affinity of the plasma protein with the compound).&#x20;

* For PPB values> 80%, it presents good affinity.
* For PPB values<40%, it presents a low level of affinity.

This is a regression model that receives SMILES as input and returns PPB values (float format) as output.

The model is pre-trained, and the author provides us with the checkpoints to be downloaded.

This model also poses a job of optimizing the Plasma Protein Binding property, therefore we are going to focus only on the part in which we obtain the **PPB values**.

### **3.  Open a model request**

To begin incorporating a model, first open a Model Request [issue](https://github.com/ersilia-os/ersilia/issues) on the Ersilia Repository.

1. Fill in all the required information such as Model Name, Model Description, Slug etc.
2. Once you submit the issue, edit it and add a new line so that the information you supplied is parsed automatically.&#x20;

### 4.  Run the model outside Ersilia.

To incorporate a model into Ersilia, we must ensure that the model provided by the third party can be executed in our local environment. To do this we must identify the place where the model code is located. For this example, we have located the [IDL-PPBopt](https://github.com/Louchaofeng/IDL-PPBopt) model in the GitHub [repository](https://github.com/Louchaofeng/IDL-PPBopt).

Most models come with a "[**README.md**](https://github.com/Louchaofeng/IDL-PPBopt#readme)" with the installation instructions, we must read them and follow the steps indicated. We must locate the code or script to run the model and the checkpoints file of the model.

To incorporate a model into the Ersilia Model Hub it's important to note the following:

* Currently Ersilia is set up to run models that require python versions 3.7 and 3.8.
* When installing the dependencies, it is recommended to do it within a conda environment created from the **CLI**. Review the basic requirements, such as ML and python packages needed to run the model and make adjustments if necessary, thus ensuring that it functions with those specific Python versions.
* Ersilia only allows running models on the **CPU**. If we find a model that uses the **Pytorch library**, keep in mind that its configuration and installation must be supported for **CPU** and not for **GPU**.

Following the example of the “[IDL-PPBopt”](https://github.com/Louchaofeng/IDL-PPBopt) model, these are the steps taken to test the model.

#### 4.1.  Identification of the list of requirements to install in the conda environment.

For this model, the official documentation from the repository lists these dependencies as requirements for executing the model.

* pytorch 1.5.0
* openbabel 2.4.1
* rdkit
* scikit learn
* scipy
* cairosvg

#### 4.2.  Create a development environment on your local machine.

In this example, the following actions were performed to create and activate the conda environment.

```bash
conda create -n idl-ppb python=3.7
conda activate idl-ppb
```

#### 4.3.  Install dependencies.

Once the environment is activated, install the necessary packages of the model as indicated by the documentation.&#x20;

{% hint style="warning" %}
Note: If the official documentation does not specify the versions of the packages to install, these versions must be specified at the time of installation as follows:
{% endhint %}

```bash
pip install rdkit==2022.9.4 
pip install scipy==1.7.3 
pip install scikit-learn==1.0.2
pip install pandas==1.3.5
pip install matplotlib==3.5.3
conda install openbabel=2.4.1 -c conda-forge
conda install pytorch==1.5.0 torchvision==0.6.0 cpuonly -c pytorch
```

In the **IDL-PPB** example model, we need to specify the Pytorch version to be CPU only. Following the steps of the official pytorch documentation [here](https://pytorch.org/get-started/locally/#supported-linux-distributions).

To know the packages installed in the local environment and the versions we can list the packages that were installed like this: &#x20;

```
pip list 
conda list
```

{% hint style="info" %}
It is recommended to save the list of the dependencies installed correctly in the environment to a txt file for future implementations.
{% endhint %}

Another method for creating a conda environment and installing dependencies is with the use of a .yml file.

* Create a .yml file with all the required dependencies. In this example, the file was named env.yml.

```yaml
name: eos22io
channels:
- conda-forge
pytorch
dependencies:
    - python=3.7
openbabel=2.4.1
pytorch=1.5.0 
torchvision=0.6.0
pip:
scikit-learn==1.0.2
pandas==1.3.5
rdkit==2022.9.4
scipy==1.7.3
matplotlib==3.5.3

```

* With the env.yml created, a conda environment containing all the dependencies listed can be generated by running the following command.&#x20;

```bash
conda env create <name-of-environment> -f env.yml
```

This method eliminates the need to install dependencies individually.

* Activate the new conda environment

```bash
conda activate <name-of-environment>
```

#### 4.4.  Clone the repository of the original model:

In this example **IDL-PPBopt**, the model is located on a Github repository. Clone it to your local environment and follow the instructions on the [repository](https://github.com/Louchaofeng/IDL-PPBopt) to run the model. From the command console, within our workspace, we can clone the repository by running the following command.

```bash
git clone https://github.com/Louchaofeng/IDL-PPBopt.git
```

Once cloned, we can open the folder containing the model code from our development environment. In this case “[Visual Studio Code](https://code.visualstudio.com/)" is used.

#### 4.5.  Identification of the main code and complementary files to run the model:

From this model, we will focus on the code that predicts the **Plasma Protein Binding(PPB)** property of a compound.

* For this example, the [**README.md**](https://github.com/Louchaofeng/IDL-PPBopt#readme) file indicates that we can run "[IDL-PPBopt.ipynb](https://github.com/Louchaofeng/IDL-PPBopt/blob/main/Code/IDL-PPBopt.ipynb)" in the jupyter notebook. In this case, we are only interested in executing it from the command console. Therefore we obtain the file IDL-PPBopt with extension .py to execute it with the python interpreter.
* In the [**saved\_models** ](https://github.com/Louchaofeng/IDL-PPBopt/tree/main/Code/saved\_models)folder we find the checkpoints.
* An example input file with the extension .csv, [**input\_compounds.csv**](https://github.com/Louchaofeng/IDL-PPBopt/blob/main/Code/input\_compounds.csv), which contains SMILES, and will be passed as input.

#### 4.6.  Inspecting the main code “IDL-PPBopt.py”.

We are going to use the only part of this script that processes the input data, loads the model and obtains the **Plasma Protein Binding (PPB)** values of a compound.

#### 4.6.1. Read the input file:

Inside IDL-PPB.py we identify the code that processes the input file. It is the following:

{% code title="IDL-PPB.py" %}
```python
raw_filename = "input_compounds.csv"                                                             
feature_filename = raw_filename.replace('.csv','.pickle')
filename = raw_filename.replace('.csv','')
prefix_filename = raw_filename.split('/')[-1].replace('.csv','')
smiles_tasks_df = pd.read_csv(raw_filename)
smilesList = smiles_tasks_df.cano_smiles.values
```
{% endcode %}

In this example, **IDL-PPB** uses an input file named “**input\_compounds.csv**”. To read the file use Pandas and the read\_csv function. Pandas read\_csv function imports a CSV file to Dataframe format.&#x20;

Gets a list of smiles from the file and processes the smiles. In this process, it converts the smiles into a canonical format using the rdkit package.

#### 4.6.2.  Calculate the molecular features.

In this process, they calculate the bond features using the **attentive fingerprint algorithm (AFP)**.

{% code title="IDL-PPB.py" %}
```python
#Calculate the molecular feature
feature_dicts = save_smiles_dicts(smilesList,filename)
remained_df= smiles_tasks_df[smiles_tasks_df["cano_smiles"].isin(feature_dicts['smiles_to_atom_mask'].keys())]
uncovered_df = smiles_tasks_df.drop(remained_df.index)
print(str(len(uncovered_df.cano_smiles))+' compounds cannot be featured')
remained_df = remained_df.reset_index(drop=True)


```
{% endcode %}

Here you can check if there are some compounds that cannot be featured. Find more information about the **Attentive fingerprint model** [here](https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b00959).

#### 4.6.3.  Load the pre-trained model.

This model uses the [PyTorch library](https://pytorch.org/docs/stable/index.html) to load the model. We can see it in the following piece of code:

{% code title="IDL-PPB.py" %}
```python
moIDL-PPB.pydel = Fingerprint(radius, T, num_atom_features, num_bond_features, fingerprint_dim, output_units_num, p_dropout)
model.cuda()
best_model = torch.load('saved_models/model_ppb_3922_Tufe_Dec_22_22-23-22_2020_'+'54'+'.pt')
```
{% endcode %}

It can be seen that it is configured to compile with **CUDA**, which means that it uses the **GPU** to obtain the predictions. Therefore it is necessary to change it to run with **CPU**.

#### 4.6.4.  Get the prediction values.

Finally, the predictions are saved in a file named "**temp.csv**".

{% code title="IDL-PPB.py" %}
```python
remain_pred_list = eval(model, remained_df)
remained_df['Predicted_values'] = remain_pred_list
remained_df
remained_df.to_csv('temp.csv')
```
{% endcode %}

After the previous inspection of the main code to execute the model. We found that this model is coded to run on the **GPU**.&#x20;

{% hint style="warning" %}
In this case, **Ersilia** (according to the [requirements](example-of-the-model-incorporation-workflow.md#4.3.-install-dependencies.)) only allows it to be executed on the **CPU**, and the PyTorch library was installed by specifying the **CPU**.&#x20;
{% endhint %}

At this point, if we run the model in this way, we would have an error like:

{% hint style="danger" %}
```
"RuntimeError: No CUDA compatible device detected"
```
{% endhint %}

Therefore it is necessary to debug the model and remove the Cuda functionality.

#### 4.7.  Debugging the IDL-PPB model to run with CPU.

Based on PyTorch's documentation, make the changes so that the code runs correctly on the **CPU**.

In the **IDL-PPBopt.py** file, everything that is instantiated with **CUDA** must be modified. In addition to adding some lines of code to specify the device.&#x20;

The steps for this example are as follows:&#x20;

**4.7.1**.  Set the type of device to ‘**CPU**’ at the beginning of the code and after imports. View [torch.device](https://pytorch.org/docs/stable/tensor\_attributes.html?highlight=torch+device#torch.device) documentation.

<pre class="language-python" data-title="IDL-PPBopt.py"><code class="lang-python"><strong>device = torch.device( 'cpu')
</strong></code></pre>

**4.7.2.**  Change the default value of the floating point tensor type.

{% code title="IDL-PPBopt.py" %}
```python
torch.set_default_tensor_type('torch.FloatTensor')
```
{% endcode %}

**4.7.3.**  Setting the model to run on the previously configured CPU device.&#x20;

Change this code instruction:

{% code title="IDL-PPBopt.py" %}
```python
model.cuda()
```
{% endcode %}

For this:

{% code title="IDL-PPBopt.py" %}
```python
model.to(device)
```
{% endcode %}

**4.7.4.**  Change this line of code:

```
best_model = torch.load(model_pretrained)
```

**4.7.5.** To load the model to run on the CPU.&#x20;

```python
best_model=torch.load(model_pretrained,map_location=torch.device('cpu'))
```

**4.7.6**.  Everything that will instantiate cuda, such as **`torch.cuda.LongTensor,`** change it to **`torch.LongTensor`**, as indicated in the [documentation ](https://pytorch.org/docs/stable/tensors.html#torch.Tensor.cpu)to run with CPU.

**4.7.7.**  For the construction and evaluation of the model, they used the attentive fingerprint algorithm ([**AFP**](https://github.com/Louchaofeng/IDL-PPBopt/tree/main/Code/AttentiveFP)). For this, the author has implemented a class called **Fingerprint**. In this class, they also use CUDA, so everything that was as **`torch.cuda.FloatTensor`**, we should change it to **`torch.FloatTensor.`**

#### 4.8.  Test the model with the example input file.

Ensure that the model is tested thoroughly to certify that it is fully functional and returns the expected output before submitting a model request to Ersilia. Execute the command to test the model using the example input file provided.

```bash
(base) carcablop@DESKTOP-1KQGKDF:~$ conda activate idl-ppb
(idl-ppb)carcablop@DESKTOP-1KQGKDF:~$ cd IDL-PPBopt/Code/
(idl-ppb)carcablop@DESKTOP-1KQGKDF:~/IDL-PPBopt/Code$python IDL-PPBopt.py
```

The model works, and we can see an output in the terminal like this:

```
number of all smiles:  4
number of successfully processed smiles:  4
feature dicts file saved as input_compounds.pickle
0 compounds cannot be featured
/home/carcablop/miniconda3/envs/idl-ppb/lib/python3.7/site-packages/torch/serialization.py:657: SourceChangeWarning: source code of class 'torch.nn.modules.linear.Linear' has changed. Saved a reverse patch to Linear.patch. Run `patch -p0 < Linear.patch` to revert your changes.
  warnings.warn(msg, SourceChangeWarning)
/home/carcablop/miniconda3/envs/idl-ppb/lib/python3.7/site-packages/torch/serialization.py:657: SourceChangeWarning: source code of class 'torch.nn.modules.container.ModuleList' has changed. Saved a reverse patch to ModuleList.patch. Run `patch -p0 < ModuleList.patch` to revert your changes.
  warnings.warn(msg, SourceChangeWarning)
```

The result produces an output file with the smiles string and the predicted value.

| <p>cano_smiles, Predicted_values<br>0,O=C(O)CC(c1ccccc1)n1ccc2cc(OCCc3ccc4c(n3)NCCC4)ccc21,0.97072583<br>1,CN(C)Cc1cncc(C(CC(=O)O)n2ccc3cc(OCCc4ccc5c(n4)NCCC5)ccc32)c1,0.8506336<br>2,CC(C)N1CN(C(c2ccccc2)c2ccccc2)n2ccc(=O)c(O)c2C1=O,0.9469087<br>3,COCCN1CN(C(c2ccccc2)c2ccccc2)n2ccc(=O)c(O)c2C1=O,0.9236307</p> |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

#### 4.9.  Run the model with the input file provided by Ersilia.

The input file that Ersilia provides to test the models is the .csv file called "[eml\_canonical.csv](https://raw.githubusercontent.com/ersilia-os/ersilia/master/notebooks/eml\_canonical.csv)". This file has a header, and three columns like this: drugs, smiles, and can\_smiles. Unlike the sample model file, which contains a column named "cano\_smiles". \[[see point 4.6](example-of-the-model-incorporation-workflow.md#4.6.1.-read-the-input-file)] Regarding the way the author codes to read the input from the file, if we directly pass the input file "eml\_canonical.csv" an “Attribute error” will be thrown, indicating that it does not have a “cano\_smiles” attribute.&#x20;

To avoid changing the code and modifying the input file, we can add the following code:

```python
with open("eml_canonical.csv", "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[1] for r in reader]

smiles_tasks_df= pd.DataFrame({'cano_smiles': smiles_list})
```

It runs correctly and the console output is as follows:

|                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <p>number of all smiles:  442<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>[10:57:27] WARNING: not removing hydrogen atom without neighbors<br>number of successfully processed smiles:  442<br>[CaH2]<br>[10:57:31] WARNING: not removing hydrogen atom without neighbors<br>[10:57:31] WARNING: not removing hydrogen atom without neighbors<br>[F-]<br>[10:57:34] WARNING: not removing hydrogen atom without neighbors<br>[10:57:34] WARNING: not removing hydrogen atom without neighbors<br>[I]<br>[10:57:37] WARNING: not removing hydrogen atom without neighbors<br>[10:57:37] WARNING: not removing hydrogen atom without neighbors<br>O<br>[Cl-].[K+]<br>[I-].[K+]<br>S<br>[10:57:49] WARNING: not removing hydrogen atom without neighbors<br>[10:57:49] WARNING: not removing hydrogen atom without neighbors<br>N.N.[Ag+].[F-]<br>[Cl-].[Na+]<br>[10:57:55] WARNING: not removing hydrogen atom without neighbors<br>[10:57:55] WARNING: not removing hydrogen atom without neighbors<br>feature dicts file saved as eml_canonical.pickle<br>9 compounds cannot be featured<br>/home/carcablop/miniconda3/envs/idl-ppb/lib/python3.7/site-packages/torch/serialization.py:657: SourceChangeWarning: source code of class 'torch.nn.modules.linear.Linear' has changed. Saved a reverse patch to Linear.patch. Run `patch -p0 &#x3C; Linear.patch` to revert your changes.<br>  warnings.warn(msg, SourceChangeWarning)<br>/home/carcablop/miniconda3/envs/idl-ppb/lib/python3.7/site-packages/torch/serialization.py:657: SourceChangeWarning: source code of class 'torch.nn.modules.container.ModuleList' has changed. Saved a reverse patch to ModuleList.patch. Run `patch -p0 &#x3C; ModuleList.patch` to revert your changes.<br>  warnings.warn(msg, SourceChangeWarning)</p> |

#### 4.10.  Verify the model output.

The above output indicates that of the 442 compounds, there are 9 compounds for which the characteristics could not be calculated. The compounds specifically were the compounds that were printed on the console and are characterized by being inorganic and salts. The author does not take into account this type of compound, therefore it does not appear in the model output file.

{% hint style="warning" %}
Note: We must always make sure that when we test a model it must make the predictions for all the compounds in the input file "eml\_canonical.csv".
{% endhint %}

In the implementation of the model with the Ersilia template, we must take into account that for these compounds the output values must be Null.

#### 4.11.  Share model output on the GitHub issue.

Finally, the model output file should be shared in the GitHub issue of the model request. (This tells the coordinators that we've successfully installed the model's dependencies, and that the model returns the expected output.).

### **5.  Model request approval.**

Once the problems that have arisen in the installation of the original model have been resolved, and we have obtained a model output, the Ersilia coordinators will be available to approve the model.

### 6.   Receive a unique model id.

Once approved, a repository would be created with a generated Ersilia unique id. In the case of this model, the id is eos22io. This is the link to the repository of the model eos22io: [https://github.com/ersilia-os/eos22io](https://github.com/ersilia-os/eos22io)

#### 6.1.  Fork the repository and clone it to your local device.

To fork the model, go onto the GitHub page of the newly created model and click on the fork in the top right and then you can clone your forked repository.

From the ubuntu terminal. We execute the following commands:

```bash
git clone https://github.com/carcablop/eos22io.git
cd eos22io/
```

### 7.  Populate the model repo with the data and codes for prediction.

Follow the folder structure in the repo and input the necessary codes. Using the [eos22io ](https://github.com/ersilia-os/eos22io)repo as an example:

| <p>eos22io<br>│   .gitattributes<br>│   .gitignore<br>│   Dockerfile<br>│   LICENSE<br>│   metadata.json<br>│   pack.py<br>│   README.md<br>├───model<br>│   │   README.md<br>│   ├───checkpoints<br>│   │   README.md<br>│   └───framework<br>│   │   input.csv<br>│   │   README.md<br>│   │   run.sh<br>│   │   __init__.py<br>│   ├───code<br>│   │   main.py<br>│   ├───test_data_dir<br>│   │   predictions.csv<br>│   └───train_data<br>└───src<br>    service.py</p> |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |

The next step is the migration of the original code to the corresponding folders of the Ersilia template.&#x20;

According to the previous architecture, take into account the following to copy the code:

#### **7.1**.  Copy code to the Ersilia template

All the code that is needed to make the predictions should be migrated to the folder _`eos22io/model/framework/code/`_.

```bash
cd ~/Desktop
cp -r IDL-PPBopt/Code/ eos22io/model/framework/code/.
```

Now we have the corresponding code to make the predictions, and additionally, we have the **IDL-PPB.py** file that was modified in the [previous steps](example-of-the-model-incorporation-workflow.md#4.7.-debugging-the-idl-ppb-model-to-run-with-cpu.) (4.7) to run the model with the **CPU** and not with **CUDA** as it was originally coded.&#x20;

To organize the code, a folder "`idl_ppb_`" was created to store the original script and the script that was modified (**IDL-PPB.py**), additionally, the text files used as examples for the model inputs were saved in that folder.

#### 7.2.  Add the model checkpoint file to the checkpoint folder.&#x20;

This can come in several flavours depending on the original repository. For example, it could be a Pickle or Joblib file (.pkl or .joblib extensions) if the modelling framework used is Scikit Learn. If Pytorch is used to build the model, it could be a .pth, .pt, or .pth.tar file. And finally, if the model was built using Tensorflow, it could be a .ckpt file.

In this example, the checkpoints are located in the “**saved\_model**” folder provided by the author. We will save them in the checkpoints folder provided by the Ersilia template.

```bash
cp IDL-PPBopt/Code/saved_models/model_ppb_3922_Tue_Dec_22_22-23-22_2020_54.pt eos22io/model/checkpoints/.

```

#### 7.4.  Editing the ‘[main.py](https://github.com/ersilia-os/eos-template/blob/main/model/framework/code/main.py)’ of the [Ersilia template](https://github.com/ersilia-os/eos-template).

As seen in the previous steps, the file we used to run the model from the command console was the **IDL-PPB.py** file, but the Ersilia template provides us with the **main.py** file, which is the main file to run the model.&#x20;

In this example, we have decided to copy from the IDL-PPB.py file the functions to load the model and calculate the predictions (since it contains all the modifications so that the model runs without errors, [see point 4.7)](example-of-the-model-incorporation-workflow.md#4.7.-debugging-the-idl-ppb-model-to-run-with-cpu.), leaving this file (**IDL-PPB.py**) only with the functions that we are going to call from other modules of the project. The rest of the functions can be called through the imports in the **main.py**.&#x20;

The **main.py** is as follows:

We will import all the modules and packages needed to run the model:&#x20;

{% code title="main.py" %}
```python
import csv
import copy
import os
import sys
import pandas as pd
from rdkit import Chem
from scipy import stats
import torch
import torch.nn as nn
torch.manual_seed(8) # for reproduce
CUDA_VISIBLE_DEVICES = 0
sys.setrecursionlimit(50000)
torch.backends.cudnn.benchmark = True
torch.set_default_tensor_type('torch.FloatTensor')
torch.nn.Module.dump_patches = True

#then import my own modules
from AttentiveFP import Fingerprint, Fingerprint_viz, save_smiles_dicts, get_smiles_array
from sarpy.SARpytools import *
from idl_ppb_.idl_ppb_modular import *
```
{% endcode %}

#### Reading command line arguments:

{% code title="main.py" %}
```python
#parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
```
{% endcode %}

We are going to establish the relative paths of the checkpoints that we copied in the previous point.

In the case of the **eos22io** model, the route was established like this:

{% code title="main.py" %}
```python
# current file directory
root = os.path.dirname(os.path.abspath(__file__))
# checkpoints directory
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))
model_pretrained=os.path.join(checkpoints_dir, "model_ppb_3922_Tue_Dec_22_22-23-22_2020_54.pt")
```
{% endcode %}

**Set global variables:**

{% code title="main.py" %}
```python
batch_size = 64
radius = 2
T = 2
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
```
{% endcode %}

We are going to read the input file that contains the smiles, which will be used as input for the function **my\_model ()**.

{% code title="main.py" %}
```python
# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model
output=my_model(smiles_list)
```
{% endcode %}

In the my\_model function, we process the input, load the model using the PyTorch library, and calculate the predictions:

```python
# my model
def my_model(smiles_list):
    df = pd.DataFrame({'can_smiles': smiles_list})
    smiles_tasks_df = df.copy()
    smilesList = smiles_tasks_df.can_smiles.values
    remained_smiles=smiles_list

    #smiles into canonical form
    canonical_smiles_list= smiles2canonical(smilesList)
       
    smiles_tasks_df = smiles_tasks_df[smiles_tasks_df["can_smiles"].isin(remained_smiles)]
    smiles_tasks_df['can_smiles'] =canonical_smiles_list
    assert canonical_smiles_list[0]==Chem.MolToSmiles(Chem.MolFromSmiles(smiles_tasks_df['can_smiles'][0]), isomericSmiles=True)

    #Calcule the molecula feature
    feature_dicts = save_smiles_dicts(smilesList,filename)
    remained_df = smiles_tasks_df[smiles_tasks_df["can_smiles"].isin(feature_dicts['smiles_to_atom_mask'].keys())]
    uncovered_df = smiles_tasks_df.drop(remained_df.index)
    print(str(len(uncovered_df.can_smiles))+' compounds cannot be featured')
    remained_df = remained_df.reset_index(drop=False)
   
    #Load the model
    p_dropout= 0.1
    fingerprint_dim = 200
    weight_decay = 5 # also known as l2_regularization_lambda
    learning_rate = 2.5
    output_units_num = 1 # for regression model
    x_atom, x_bonds, x_atom_index, x_bond_index, x_mask, smiles_to_rdkit_list = get_smiles_array([canonical_smiles_list[0]],feature_dicts)
    num_atom_features = x_atom.shape[-1]
    num_bond_features = x_bonds.shape[-1]
    loss_function = nn.MSELoss()
    model = Fingerprint(radius, T, num_atom_features, num_bond_features,
                fingerprint_dim, output_units_num, p_dropout)
   
    model.to(device) ####change to CPU
    #model.cuda()
    best_model = torch.load(model_pretrained, map_location=torch.device('cpu')) ###change adding map_location
    #best_model = torch.load(model_pretrained)
    best_model_dict = best_model.state_dict()
    best_model_wts = copy.deepcopy(best_model_dict)
    model.load_state_dict(best_model_wts)
    (best_model.align[0].weight == model.align[0].weight).all()
    model_for_viz = Fingerprint_viz(radius, T, num_atom_features, num_bond_features,
                fingerprint_dim, output_units_num, p_dropout)
    #model_for_viz.cuda()
    model_for_viz.to(device) ####change to cpu
    model_for_viz.load_state_dict(best_model_wts)
    (best_model.align[0].weight == model_for_viz.align[0].weight).all()

    #Predict values
    remain_pred_list = eval(model, remained_df,feature_dicts)
    remained_df['Predicted_values'] = remain_pred_list
```

Taking into account what was mentioned in [previous steps](example-of-the-model-incorporation-workflow.md#4.10.-verify-the-model-output.) ([step 4.10](example-of-the-model-incorporation-workflow.md#4.10.-verify-the-model-output.)), for some compounds such as salts and inorganics that do not predict the PPB values, it was decided that for these compounds the output should be Null. Therefore, the following function is implemented:

```python
#making sure it returns all values.For the compounds that it is not possible to calculate the features the output is null
    array_ppbs= remained_df['Predicted_values'].values
    array_missindex= uncovered_df.index.values
    print(array_missindex)
    for indice in array_missindex:
        if indice < len(array_ppbs):
            array_ppbs = np.insert(array_ppbs, indice, None)
    return array_ppbs
```

Finally we return the list of PPB values that it predicts for all the compounds we pass as input to the model.

The model output is written in a .csv format

```python
# write PPB values output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow(["PPB_VALUES"])  # header
    for o in output:
        writer.writerow([o])
```

#### 7.5.  Edit the service.py:&#x20;

The service.py file is only modified if necessary.&#x20;

Only this part of the code: [service.py](https://github.com/ersilia-os/eos-template/blob/main/src/service.py#L93) will be modified to specify the type of output, whether it is String or Float type.

In this case, the PPB values are of Float type. And in the service.py file, the output type is already set as Float, therefore it was not necessary to modify it.

#### 7.6.  Edit the Dockerfile&#x20;

With all the dependencies required to run the model. This is essential for the successful incorporation of the model into Ersilia.

We will copy the same dependencies used in [step 4.3](example-of-the-model-incorporation-workflow.md#4.3.-install-dependencies.):

```docker
FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN pip install rdkit==2022.9.4
RUN pip install scipy==1.7.3
RUN pip install scikit-learn==1.0.2
RUN pip install pandas==1.3.5
RUN pip install matplotlib==3.5.3
RUN conda install openbabel=2.4.1 -c conda-forge
RUN conda install pytorch==1.5.0 torchvision==0.6.0 cpuonly -c pytorch

WORKDIR /repo
COPY . /repo
```

Since model files can be large and it is generally not a good practice to track large files within a Git repo, model checkpoints are tracked on a **Git LFS** (Large File Storage) server. Make sure you have Git LFS installed and initialized on your development environment. Run the following command to track the model checkpoints:

```git
git lfs track model/checkpoints/<model-checkpoint>
```

This will create a new entry in the **.gitattributes file**. Commit this change and push it to the repository. Now you can track everything using Git normally and continue on the next steps.

#### 7.7.  Edit the [metadata.json](https://github.com/ersilia-os/eos-template/blob/main/metadata.json) file

Following the instructions specified in the Anatomy of the Ersilia Model Template section of the Ersilia [gitbook](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/contribute-models/model-incorporation-guidelines). Make sure to match the case of accepted values exactly. For example, if the relevant tags for your model are “ADME”, you should write “ADME” and not “adme”, or “Adme”. Similarly if the , tag is “Bioactivity profile”, you should write it as such and not as “bioactivity profile”, or “Bioactivity Profile”.

```json
{    
    "Identifier": "eos22io",
    "Slug": "idl-ppbopt",
    "Status": "In progress",
    "Title": "Human Plasma Protein Binding (PPB) of Compounds",
    "Description": "IDL-PPB aims to obtain the plasma protein binding (PPB) values of a compound. Based on an interpretable deep learning model and using the algorithm fingerprinting (AFP) this model predicts the binding affinity of the plasma protein with the compound.",
    "Mode": " Pretrained",
    "Task": ["Regression "],
    "Input": [" Compound"],
    "Input Shape": "Single ",
    "Output": ["Experimental value"],
    "Output Type": [" Float"],
    "Output Shape": "Single ",
    "Interpretation": "This model receives smiles as input and returns as output the fraction PPB, which measures the affinity of the binding of the plasma protein. In the analysis of results by the author, they indicate high affinity (fraction of ppb >80%), medium affinity (40% <= fraction of ppb <=80%) and as low levels of afity (fraction of ppb < 40%). Note: Inorganics and salts are out of the applicability domain of the model, So for these compounds the output is Null. ",
    "Tag": [
        "Fraction bound","ADME"
    ],
    "Publication": "https://pubs.acs.org/doi/10.1021/acs.jcim.2c00297",
    "Source Code": "https://github.com/Louchaofeng/IDL-PPBopt",
    "License": "GPL-3.0"
}
```

If the [accepted tags](https://github.com/ersilia-os/ersilia/blob/master/ersilia/hub/content/metadata/tag.txt) do not contain a tag relevant to your model, open a PR with this tag added to the file ersilia/hub/content/metadata/tags.txt within the [ersilia repository](https://github.com/ersilia-os/ersilia).&#x20;

Edit the **README** file in the `model/framework/code` directory to outline the installation instructions to set up the development environment locally.

#### 7.8.  Run the main.py file and generate the prediction results.

Use a sample input file to generate an example output file with the data from the model.

First, **we must activate the conda environment,** where all the installed model dependencies are located ([see step 4.2](example-of-the-model-incorporation-workflow.md#4.2.-create-a-development-environment-on-your-local-machine.)).

```python
cd model/framework/code
python main.py ../<name_of_input_file>.csv ../<name_of_output_file>.csv
```

The model runs correctly and we get the output file without errors. The output corresponds to the expected value

### 8.  Run the model within the Ersilia CLI.

The next step involves testing the model within the Ersilia CLI.&#x20;

To execute this step, Ersilia must be correctly installed on your system, see the steps here:  [gitbook](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/installation).

Activate the ersilia environment.

In the command console, we execute the following command.

```bash
conda activate ersilia
```

Once the ersilia environment is active, we can test the model within the Ersilia CLI using the path to your clone of the repo on your local device i.e. `c/Users/<username>/Desktop/eos22io`

```bash
ersilia -v fetch model_id -r <local_path_to_repo>
```

In the command console we can see it like this:

```bash
(ersilia) user_name@DESKTOP-1KQGKDF:~$ ersilia -v fetch eos22io --repo_path /home/user_name/eos22io/
```

If there were no issues encountered during this step, a success message would be returned. This shows that the model is ready to be incorporated into Ersilia.

```
👍 Model eos22io fetched successfully!
```

**Serve the model:**

```
ersilia -v serve eos22io
```

**Make predictions**

{% code title="" %}
```
ersilia -v api run -i "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]"
```
{% endcode %}

{% code title="# save output in a CSV file" %}
```
ersilia -v api run -i eml_canonical.csv -o output_eos22io.csv
```
{% endcode %}

#### 8.1.  Share model output on the GitHub issue.

Finally the model output file should be shared in the GitHub issue of the model on GitHub.&#x20;

#### 8.2.  Update the changes to your fork of the model repository.

To publish the changes on the cloned repository, all the required files need to be added to the GitHub repository. Use this command to add all the model files to git, create a commit message with the information for version control and push the code to the main branch of the remote repository.

```
git commit -am "initial commit" && git push origin main
```

#### 8.3.  Make a pull request to merge changes to the Ersilia repository.

The next step in the model incorporation is creating a pull request to merge the changes in your repository with the Ersilia repository.

After pushing your codes to your remote model repository, navigate to the main branch of the repository on GitHub and make contributions to making a pull request to the main branch of the Ersilia model repository

#### 8.4.  The model must be tested by others.

The model must be successfully tested by other collaborators both in the CLI and in Google Colab, and the logs and the output files of the model must be shared. Finally, the model can be incorporated into the Ersilia model Hub.

### 9. Deleted Fork of the model

Delete the fork from your personal repository.
