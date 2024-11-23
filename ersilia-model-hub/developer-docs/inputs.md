---
description: Here we describe the types of valid inputs while running Ersilia models
---

# Inputs

## Compound inputs

Small molecule structures are typically expressed in the well known [SMILES](https://en.wikipedia.org/wiki/Simplified\_molecular-input\_line-entry\_system) format. Below are the SMILES strings of a few drug molecules:

<table><thead><tr><th width="166">Drug</th><th>SMILES</th></tr></thead><tbody><tr><td>Artemisin</td><td>CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O</td></tr><tr><td>Isoniazid</td><td>C1=CN=CC=C1C(=O)NN</td></tr><tr><td>Tenofovir</td><td>CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O</td></tr><tr><td>Aspirin</td><td>CC(=O)OC1=CC=CC=C1C(=O)O</td></tr><tr><td>Ibuprofen</td><td>CC(C)CC1=CC=C(C=C1)C(C)C(=O)O</td></tr><tr><td>Remdesivir</td><td>CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C</td></tr><tr><td>Cephalotaxin</td><td>COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5</td></tr></tbody></table>

In the most common case, each input sample corresponds to a **single molecule**. However, some models expect a **list of molecules** as input, and some even expect more complex inputs such as **pairs of lists of molecules**. Multiple molecules (lists) can be serialized to strings in the SMILES notation with the dot (`.`) character.

Valid **input file formats** are comma-separated (`.csv`), tab-separated (`.tsv`)  and JSON (`.json`). Ersilia automatically recognizes these formats. Inputs can also be passed as Python instances through the **Ersilia Python API**.

It is possible to run Ersilia models for **one input** as well as **multiple inputs**. Ersilia automatically detects the if one or multiple inputs are passed.

### Single molecules

#### One input

This is the simplest case where one single molecule is passed as input.

{% tabs %}
{% tab title="CSV" %}
{% code title="compound_single.csv" %}
```csv
smiles
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O
```
{% endcode %}
{% endtab %}

{% tab title="JSON" %}
{% code title="compound_single.json" %}
```json
"CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
```
{% endcode %}
{% endtab %}

{% tab title="Python" %}
```python
smiles = "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
```
{% endtab %}
{% endtabs %}

#### Multiple inputs

This is a common case too, where multiple single molecules are passed as input. The model will run predictions/calculations for each molecule independently.

{% tabs %}
{% tab title="CSV" %}
{% code title="compound_singles.csv" %}
```csv
smiles
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O
C1=CN=CC=C1C(=O)NN
CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O
CC(=O)OC1=CC=CC=C1C(=O)O
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C
COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5
```
{% endcode %}
{% endtab %}

{% tab title="JSON" %}
{% code title="compound_singles.json" %}
```json
[
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
    "C1=CN=CC=C1C(=O)NN",
    "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
    "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5"
]
```
{% endcode %}
{% endtab %}

{% tab title="Python" %}
```python
smiles = [
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
    "C1=CN=CC=C1C(=O)NN",
    "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
    "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5",
]
```
{% endtab %}
{% endtabs %}

{% hint style="info" %}
The majority of Ersilia chemistry models take single molecules as input.
{% endhint %}

### List of molecules

#### One input

Here the molecule expects a list of molecules as input, therefore, the one prediction/calculation will be done based on the list as a whole. Some generative models, for example, require multiple molecules as a starting point for one generation round.

{% tabs %}
{% tab title="CSV" %}
{% code title="compound_list.csv" %}
```
smiles
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O
C1=CN=CC=C1C(=O)NN
CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O
CC(=O)OC1=CC=CC=C1C(=O)O
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C
COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5
```
{% endcode %}
{% endtab %}

{% tab title="JSON" %}
{% code title="compound_list.json" %}
```json
[
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
    "C1=CN=CC=C1C(=O)NN",
    "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
    "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5"
]
```
{% endcode %}
{% endtab %}

{% tab title="Python" %}
```python
smiles_list = [
    "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
    "C1=CN=CC=C1C(=O)NN",
    "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    "CC(=O)OC1=CC=CC=C1C(=O)O",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
    "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5",
]
```
{% endtab %}
{% endtabs %}

{% hint style="warning" %}
Note that, in this case, the `compound_list.csv` file is the same as the `compound_singles.csv` file, corresponding to the multiple single molecules. However, the model will treat these files differently. In the current case, the full list corresponds to one input, whereas in the previous case each single molecule was an independent input. In the [Ersilia Model Hub](https://airtable.com/shrUcrUnd7jB9ChZV/tblZGe2a2XeBxrEHP), the **Input Shape** field is labelled as _Single_ or _List_, correspondingly.
{% endhint %}

#### Multiple inputs

Here, multiple lists are passed as input, and each list is treated independently by the model.

{% tabs %}
{% tab title="CSV" %}
{% code title="compound_lists.csv" %}
```csv
smiles
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O.C1=CN=CC=C1C(=O)NN.CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O
CC(=O)OC1=CC=CC=C1C(=O)O.CC(C)CC1=CC=C(C=C1)C(C)C(=O)O.CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C.COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5
```
{% endcode %}
{% endtab %}

{% tab title="JSON" %}
{% code title="compound_lists.json" %}
```json
[
    [
        "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
        "C1=CN=CC=C1C(=O)NN",
        "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O"
    ],
    [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
        "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5"
    ]
]
```
{% endcode %}
{% endtab %}

{% tab title="Python" %}
```python
smiles_lists = [
    [
        "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
        "C1=CN=CC=C1C(=O)NN",
        "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    ],
    [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
        "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5",
    ],
]
```
{% endtab %}
{% endtabs %}

{% hint style="info" %}
To specify a list in a single column in a `.csv` file, molecules can be separated with a dot (`.`). **The type of delimiter is specific to the input type**. The SMILES notation naturally accepts the dot as a separator for multiple molecules.
{% endhint %}

### Pair of lists of molecules

#### One input

This corresponds to a less common case where, one input is expressed as a pair of lists. An example would be a model comparing two sets of molecules and returning an overall similarity values (one float number) between the two sets.

{% tabs %}
{% tab title="CSV" %}
{% code title="compound_pair_of_lists.csv" %}
```csv
smiles_1,smiles_2
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O,CC(=O)OC1=CC=CC=C1C(=O)O
C1=CN=CC=C1C(=O)NN,CC(C)CC1=CC=C(C=C1)C(C)C(=O)O
CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O,CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C
,COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5
```
{% endcode %}
{% endtab %}

{% tab title="JSON" %}
{% code title="compound_pair_of_lists.json" %}
```json
[
    [
        "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
        "C1=CN=CC=C1C(=O)NN",
        "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O"
    ],
    [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
        "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5"
    ]
]
```
{% endcode %}
{% endtab %}

{% tab title="Python" %}
```python
smiles_pair_of_lists = (
    [
        "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
        "C1=CN=CC=C1C(=O)NN",
        "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    ],
    [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
        "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5",
    ],
)
```
{% endtab %}
{% endtabs %}

{% hint style="warning" %}
Please note that the `compound_pair_of_lists.csv` file contains two columns, one for each set. The first set has three molecules, and the second set has four molecules. Therefore, the first column contains one empty row.
{% endhint %}

#### Multiple inputs

Multiple pairs of lists can be passed to obtain multiple predictions/calculations, one for each pair. Like in the case of multiple lists, molecules can be separated with a dot character in a tabular file.

{% tabs %}
{% tab title="CSV" %}
{% code title="compound_pair_of_lists.csv" %}
```csv
smiles_1,smiles_2
CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O.C1=CN=CC=C1C(=O)NN,CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O
CC(=O)OC1=CC=CC=C1C(=O)O.CC(C)CC1=CC=C(C=C1)C(C)C(=O)O,CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C.COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5
```
{% endcode %}
{% endtab %}

{% tab title="JSON" %}
{% code title="compound_pairs_of_lists.json" %}
```json
[
    [
        [
            "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O",
            "C1=CN=CC=C1C(=O)NN"
        ],
        [
            "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O"
        ]
    ],
    [
        [
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
        ],
        [
            "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
            "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5"
        ]
    ]
]
```
{% endcode %}
{% endtab %}

{% tab title="Python" %}
```python
smiles_pairs_of_lists = [
    [
        ["CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O", "C1=CN=CC=C1C(=O)NN"],
        ["CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O"],
    ],
    [
        ["CC(=O)OC1=CC=CC=C1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"],
        [
            "CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C",
            "COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5",
        ],
    ],
]
```
{% endtab %}
{% endtabs %}

### Real-world sets of molecules

You can generate inputs of arbitrary size for your model of interest with the following command. In this case, we generate 1,000 inputs for the `chemprop-antibiotic` model and store them as a `.csv` file.

```bash
ersilia example chemprop-antibiotic -n 1000 -f input.csv
```

This command simply samples drug molecules from the attached table.

{% file src="../.gitbook/assets/drug_molecules.tsv" %}
Tab-separated file containing drug molecules from the [Drug Repurposing Hub](https://www.broadinstitute.org/drug-repurposing-hub). InChIKeys, SMILES and names are provided.
{% endfile %}

{% hint style="info" %}
Ersilia will automatically detect the SMILES column and the format in an input file, so it is acceptable to pass the `drug_molecules.tsv` file as is, or a chunk of it.
{% endhint %}

#### Small molecule databases

Many small molecule databases exist in the public domain. If you want to look for a molecule of interest, consider the following resources:

* [PubChem](https://pubchem.ncbi.nlm.nih.gov/) as a go-to search tool to obtain generalistic chemical information.
* [ChEMBL](https://www.ebi.ac.uk/chembl/) as a search tool for bioactivity data of medicinal chemistry compounds.
* [DrugBank](https://go.drugbank.com/) to obtain comprehensive information about drug molecules.
* [ZINC](https://zinc20.docking.org/) to search for commercially-available libraries of compounds.
