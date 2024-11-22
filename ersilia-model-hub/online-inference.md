---
description: Documentation to run Ersilia models online
---

# Online inference

To ensure users from all backgrounds are able to benefit from our tools, we provide a ready-to-use, no-code solution to obtain predictions for your molecules.&#x20;

## 1. Select your model of interest

We offer a broad range of models, from bioactivity prediction against several pathogens (malaria, tuberculosis, schistosomiasis, ESKAPE pathogens...) to ADME endpoints and toxicity predictions. Use our [dynamic interface](https://ersilia.io/model-hub) to browse models according to your needs and take note of the model identifier you wish to use!

{% hint style="success" %}
Please ensure you understand the output of each model before using it. _For example, a classifier will output the probability of Active for a particular assay, i.e probability that a molecule kills the malaria parasite in an in vitro assay, and a regressor might output the predicted IC50 value at which the molecule inhibits the growth of the parasite._&#x20;
{% endhint %}

## 2. Prepare your input data

The molecules must be displayed in SMILES notation. You can use [PubChem](https://pubchem.ncbi.nlm.nih.gov/) to find the SMILES notation of a given compound, simply look for the compound name on the search bar (for example, aspirin), select the best result and scroll to the SMILES section under Name and Identifier (in this case; CC(=O)OC1=CC=CC=C1C(=O)O). If your starting input data is an `.sdf` file, use your preferred visualiser, like ChemDraw, to open the molecule and obtain its SMILES representation. You can also use free software like [Marvin.js](https://marvinjs-demo.chemaxon.com/latest/demo.html) to draw a molecule and then simply click on save 💾 it as a SMILES.

Collect your list of SMILES in a `.csv` or `.txt` file.

## 3. Run predictions & download results

Go to our [online inference app](https://ersilia-self-service.streamlit.app/) and select your model of choice from the drop down list. Copy the list of SMILES (maximum allowed 100 molecules) and click on Run Predictions. Wait a few minutes to download your results!

{% hint style="info" %}
If you wish to run larger annotations, for example running several predictions against a database of >1k molecules, please contact Ersilia directly to obtain a customised solution: hello@ersilia.io
{% endhint %}

## 4. Check your predictions

By default, Ersilia will return a `.csv` file containing the following columns:

* SMILES: the input SMILES (note that these might have been standardised if they were not provided in the standard format)
* InChIKey: 27-character unique identifier of a molecule based on the International Chemical Identifier
* Model output: one or several columns containing the predictions of the selected model. Make sure to read about the model in the literature or Ersilia documentation to appropriately interpret the model's results.

{% hint style="danger" %}
Posting to this free online service will make your molecules public. Please consider [local inference](local-inference.md) if you are working with IP-protected molecules.
{% endhint %}
