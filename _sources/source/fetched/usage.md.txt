# Quick start

**This section is work in progress! Please don't try to run code as of yet. We should be done in the following weeks.**

## Select your model

In this quick-start guide, we will run antibiotic activity predictions for small molecule compounds as an example.

You can browse available models from the [Ersilia model hub](https://ersilia-os.github.io/ersilia-hub.github.io/). Browse by simply using keyword search or define the model you are interested in by using the following categories
- input: type of data you have _molecule, DNA, RNA, protein, image, text_
- endpoint: prediction of your interest _toxicity, drug-target interaction, cell sensitivity, synthetic accessibiity..._
- application: focus on a certain disease, if applicable _cancer, malaria, HIV, COVID19..._

In this example, the model is classified as:
- input: molecule
- endpoint: antibiotic activity
- application: infectious disease

Models in Ersilia are identified with the following ID structure: `eos{digit}{3 alphanumeric characters}`. Please, use this ID when citing our models.

In this example, ID for the [antibiotic activity](https://ersilia-os.github.io/ersilia-hub.github.io/first-gemma-post/) model is eos0aaa

Our model hub includes both models developed by Ersilia as well as models available from the literature (third party authorship). In the latter, Ersilia re trains the models if needed and bundles them for easy deployment. Appropriate author citations are always included.

The antibiotic predictor is based on [Stokes et al. Cell, 2020](https://pubmed.ncbi.nlm.nih.gov/32084340/). Authors kindly shared to the community the model checkpoint.

## Run on a web app

Once you identify the models you want to use, simply click on them to find a summary of the most important features (e.g how many datapoints were used for training, whether it was experimentally validated, its application domain...) It is important to keep in mind the Limitations of each model. At Ersilia we do our best to summarize the most important, but it is up to each user to understand their data and the scope of the chosen models. If you need help, please contact us at [hello@ersilia.io](mailto:hello@ersilia.io).

For eos0aaa we have created the following summary and limitations list:

Summary
- Predicts **antibiotic activity**
- Takes **compound structures** as input
- Trained with **experimental** bioactivity data against E.coli
- Based on a dataset of **>2,000** experiments
- Results validated **experimentally**
- Used for **drug repurposing**
- Identified a **novel broad-spectrum** antibiotic
- Published in _Stokes et al., Cell 2020_

Limitations
- Training set using only _E.coli_ data, not other organisms
- Predicts antibiotic activity against a pathogen, not valid for host targets
- Ersilia retrained the model on-premise using the published datasets

To use our models, you can query a single input at a time using our website app or, if you are familiar with coding, you can download and install it in your computer following the installation instructions.

To open the website app, go to the [button] at the top of the model card.

## Run programmatically

### Command-line interface (CLI)

Ersilia models are stored as BentoML bundles. Hence one can use the outstanding functionalities of this library to deploy models.

To fetch a BentoML bundle from our repository using the CLI, simply run:

```bash
ersilia fetch eos0aaa
```

Now, you can serve predictions the way you wish. BentoML offers a great spectrum of possibilities and Ersilia incorporates all of them.
For example, you can simply run to obtain predictions for the drug Halicin:

```bash
ersilia predict eos0aaa --input='["C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]"]'
```

### Python console

#### Fetch and initialize model

```python
from ersilia import ErsiliaModel
model_id = 'eos0aaa'
em = ErsiliaModel(model_id)
```

What happens behind this code is the following:

1. The model of interest (`eos0aaa`) is fetched from:
    - A [GitHub repository](https://github.com/ersilia-os/eos0aaa) containing, _at least_:
        - A `service.py` script defining the serving class (i.e. the class that will provide predictions)
        - A `pack.py` script that will generate a model bundle for easy deployment, based on [BentoML](https://www.bentoml.ai/)
    - A compressed file, stored at [Open Science Framework (OSF)](https://osf.io/hu3km/) containing model checkpoint files, i.e. the trained parameters of the model
2. A `~/bentoml/repository/eos0aaa/` folder is created containing the bundled model (see BentoML docs for more details)
3. A `pip` package named `eos0aaa` is created in your local computer

The fetching functionalities of Ersilia are configurable. Please see advanced usage cases for more details.

#### Run prediction

Once the model is ready, one can simply run predictions for the molecules of interest. Let's say that we want to predict the antibiotic activity of:
The standard input are SMILES strings.

- Ceftazidime: `[O-]C(=O)C1=C(CS[C@]2([H])[C@H](NC(=O)C(=N/OC(C)(C)C(O)=O)\C3=CSC(N)=N3)C(=O)N12)C[N+]1=CC=CC=C1`
- Halicin: `C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]`
- Ibuprofen: `CC(C)CC1=CC=C(C=C1)C(C)C(O)=O`

```python
smiles = ['[O-]C(=O)C1=C(CS[C@]2([H])[C@H](NC(=O)C(=N/OC(C)(C)C(O)=O)\C3=CSC(N)=N3)C(=O)N12)C[N+]1=CC=CC=C1',
          'C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]',
          'CC(C)CC1=CC=C(C=C1)C(C)C(O)=O']
em.predict(smiles)
```
