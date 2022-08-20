---
description: Here we describe the types of valid inputs while running Ersilia models
---

# Inputs

## Compound inputs

Small molecule structures are typically expressed in the well known [SMILES](https://en.wikipedia.org/wiki/Simplified\_molecular-input\_line-entry\_system) format. Below are the SMILES strings of a few drug molecules:

| Drug         | SMILES                                        |
| ------------ | --------------------------------------------- |
| Artemisin    | CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O         |
| Isoniazid    | C1=CN=CC=C1C(=O)NN                            |
| Tenofovir    | CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O          |
| Aspirin      | CC(=O)OC1=CC=CC=C1C(=O)O                      |
| Ibuprofen    | CC(C)CC1=CC=C(C=C1)C(C)C(=O)O                 |
| Remdesivir   | CC1(OC2C(OC(C2O1)(C#N)C3=CC=C4N3N=CN=C4N)CO)C |
| Cephalotaxin | COC1=CC23CCCN2CCC4=CC5=C(C=C4C3C1O)OCO5       |

In the most common case, each input sample corresponds to a **single molecule**. However, some models expect a **list of molecules** as input, and some even expect more complex input such as **pairs of lists of molecules**. Multiple molecules (lists) can be serialized to strings in the SMILES notation with the dot (`.`) character.

Valid **file formats** are comma-separated (`.csv`), tab-separated (`.tsv`)  and JSON (`.json`). Inputs can also be passed as Python instances through the **Python API**.

It is possible to run Ersilia models for **one input** as well as **multiple inputs**. Ersilia automatically detects the type of input.

Below we show a few exemplary&#x20;

### Single molecules

The majority of chemistry models take as input a single molecule.

#### One input

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
{% code title="compound_single.py" %}
```python
smiles = "CC1C2C(CC3(C=CC(=O)C(=C3C2OC1=O)C)C)O"
```
{% endcode %}
{% endtab %}
{% endtabs %}

#### Multiple inputs

### List of molecules

#### One input



#### Multiple inputs



### Pair of lists of molecules



### Generate your own inputs

You can generate inputs of arbitrary size for your model of interest with the following command. In this case, we generate 1,000 inputs for the `chemprop-antibiotic` model and store them as a `.csv` file.

```bash
ersilia example chemprop-antibioitic -n 1000 -f input.csv
```

## Outputs
