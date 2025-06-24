---
description: Documentation to run Ersilia models online
---

# Online inference

To ensure users from all backgrounds are able to benefit from our tools, we provide a ready-to-use, no-code solution to obtain predictions for your molecules.&#x20;

## 1. Select your model of interest

We offer a broad range of models, from bioactivity prediction against several pathogens (malaria, tuberculosis, schistosomiasis, ESKAPE pathogens...) to ADME endpoints and toxicity predictions. Use our [dynamic interface](https://ersilia.io/model-hub) to browse models according to your needs and take note of the model identifier you wish to use!

{% hint style="success" %}
Please ensure you understand the output of each model before using it. _For example, a classifier will output the probability of being active in a particular assay, e.g. probability that a molecule kills the malaria parasite in an in vitro assay, and a regressor might output the predicted IC50 value at which the molecule inhibits the growth of the parasite._&#x20;
{% endhint %}

## 2. Prepare your input data

The molecules must be displayed in SMILES notation. You can use [PubChem](https://pubchem.ncbi.nlm.nih.gov/) to find the SMILES notation of a given compound: simply introduce the compound name on the search bar (for example, aspirin), select the best result and scroll to the SMILES section within "Name and Identifiers" (in this case; CC(=O)OC1=CC=CC=C1C(=O)O). If your starting input data is an `.sdf` file, use your preferred visualiser, like ChemDraw, to open the molecule and obtain its SMILES representation. To deal with multiple molecular file formats, including SDF, you can use [OpenBabel](https://openbabel.org/) to convert them into SMILES notation. Alternatively, you can also use free software like [Marvin.js](https://marvinjs-demo.chemaxon.com/latest/demo.html) to draw a molecule and then simply click on save 💾 it as a SMILES.

Collect your list of SMILES in a `.csv` or `.txt` file **containing a single column with a header**, for example "smiles" or "input"**.**

## 3. Run predictions and download results

Go to our [online inference app](https://ersilia-self-service.streamlit.app/) and select your model of choice from the drop down list. Copy the list of SMILES (maximum allowed 100 molecules) and click on "Run Predictions!". Wait a few minutes to download your results.

{% hint style="danger" %}
Posting to this free online service will make your molecules public. Please consider [local inference](local-inference.md) if you are working with IP-protected molecules.
{% endhint %}

{% hint style="info" %}
If you wish to run larger annotations, for example running several predictions against a database of >1k molecules, please contact Ersilia directly to obtain a customised solution: hello\[at]ersilia.io
{% endhint %}

## 4. Check your predictions

By default, Ersilia will provide a downloadable `.csv` file summarizing the results, containing the following columns:

* key: 32-character unique identifier of the molecule created by Ersilia
* input: the input SMILES (please note that these might have been standardised if they were not provided in the standard format).
* Model output: one or several columns containing the predictions of the selected model. Explantion about each column's meaning can be found on the model's GitHub repository under `model/framework/columns`. For example, for model `eos3b5e` you can go to [https://github/ersilia-os/eos3b5e](https://github.com/ersilia-os/eos3b5e) and navigate to the [run\_columns.csv](https://github.com/ersilia-os/eos3b5e/blob/main/model/framework/columns/run_columns.csv) file.
