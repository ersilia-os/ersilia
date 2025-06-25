---
description: 4-day workshop on the use of AI for Drug Discovery
---

# AI2050 AI for Drug Discovery

## Course overview

Each day is organised in four main activities:

* Keynote Lectures: a presentation by H3D and Ersilia mentors and invited speakers, introducing key concepts that we will put into practice during the day
* Skills Workshop: a hands-on training where the course facilitators will walk the participants through a specific ML exercise
* Breakout: participants will be divided in small groups to complete an exercise related to the skills workshop
* Group Discussion: each group will give a short (10 min) presentation on their breakout activity.

### Course Contents <a href="#course-contents" id="course-contents"></a>

<table><thead><tr><th width="148">Activity</th><th>Session 1</th><th>Session 2</th><th>Session 3</th><th>Session 4</th></tr></thead><tbody><tr><td>Keynote</td><td><p><a href="https://drive.google.com/file/d/1AXiCDF5bbsBNYXkgCS4D4f2iyVnGZAiS/view?usp=sharing">Slides A</a></p><p><a href="https://drive.google.com/file/d/14HlsNayCnghifod9JyfXaLGBVcRlEtW8/view?usp=sharing">Slides B</a><br><a href="https://www.menti.com/alsc5obgc7z6">Menti</a></p></td><td><a href="https://drive.google.com/file/d/1bu1Q2lqNemhFbiSweMa-PYx6aWsQfa42/view?usp=drive_link">Slides A</a><br><a href="https://drive.google.com/file/d/1pSrGRyvrpmPndBGaNaQutKU31AIEs52t/view">Slides B</a></td><td><p><a href="https://drive.google.com/file/d/1fnrxPY_uY-Zz5PEWVyUwwmP0P8JnKRuS/view?usp=drive_link">Slides A</a></p><p><a href="https://drive.google.com/file/d/1WlBiqg3UCEMZ4e6bWGY1aAGG_wtt4lH0/view">Slides B</a></p></td><td><a href="https://drive.google.com/file/d/1ThL6DJmKH12bJ7RYLsQO-_wwDm3uDrBV/view?usp=drive_link">Slides A</a><br><a href="https://drive.google.com/file/d/1EW-xorln-lWQ-1vUzePMTfU1c0LBQ6qQ/view?usp=drive_link">Slides B</a></td></tr><tr><td>Skills development</td><td><a href="https://ai2050-day1.streamlit.app/">App</a><br><a href="https://drive.google.com/file/d/18uKYzlF69mUaZjEUAMN_JD7jhDal4Jqa/view?usp=drive_link">Slides</a><br><a href="https://www.menti.com/alfgx5p5x3qq">Menti</a></td><td><p><a href="https://ai2050-day2.streamlit.app/">App</a></p><p><a href="https://drive.google.com/file/d/1ewhaWpRMgru2PCNeideP3hiP1pD2dL-h/view?usp=drive_link">Slides</a><br><a href="https://www.menti.com/altn98ny9d6h">Menti</a></p></td><td><p><a href="https://ai2050-day3.streamlit.app/">App</a></p><p><a href="https://drive.google.com/file/d/1unfnDAfFrtoCpVS9TZuVvP_dulCxmSdo/view?usp=drive_link">Slides</a><br><a href="https://www.menti.com/alypgk7nwtd5">Menti</a></p></td><td><a href="https://drive.google.com/file/d/1Qwx3mSA_DdUQ7pbecsd9UgSjTdMguJYF/view?usp=drive_link">Slides</a></td></tr><tr><td>Breakout</td><td><p><a href="https://ai2050-day1.streamlit.app/">App</a></p><p><a href="https://drive.google.com/file/d/1zbgxw40dtRVrtUWyCZACtH88R8UtEECI/view?usp=drive_link">Slides</a></p><p><a href="ai2050-ai-for-drug-discovery.md#breakout-session-day-1">Guide</a></p></td><td><p><a href="https://ersilia.io/model-hub">App</a></p><p><a href="https://drive.google.com/file/d/13kyE9HKACb2hShcMFT3nYEUKyDRQ8Hf-/view?usp=drive_link">Slides</a></p><p><a href="ai2050-ai-for-drug-discovery.md#breakout-session-day-2">Guide</a></p></td><td><p><a href="https://ai2050-day3.streamlit.app/">App</a></p><p><a href="https://drive.google.com/file/d/1WlBiqg3UCEMZ4e6bWGY1aAGG_wtt4lH0/view?usp=drive_link">Slides</a><br><a href="ai2050-ai-for-drug-discovery.md#breakout-session-day-3">Guide</a></p></td><td></td></tr></tbody></table>

{% hint style="warning" %}
After the course, the apps will remain active on StreamlitCloud and also available for local deployment from GitHub. Please keep in mind that these are apps built for the course, and might not be suitable to apply directly to a research project. Contact Ersilia for more!
{% endhint %}

## Important links

In case you need further information to follow the course, be sure to check the following:

* Ersilia Model Hub: list of all [available models](https://ersilia.io/model-hub).
* Ersilia Model Hub: [self-service](https://github.com/ersilia-os/ersilia-self-service) on GitHub (requires a GitHub account).
* Ersilia Model Hub: [local installation](broken-reference) documentation (only expert users).
* [Publications](https://drive.google.com/drive/folders/1lSan_efMH-mOVz41HM-aNvAqyIHcQ0Fr?usp=sharing) referred to during the workshop.
* [Datasets](https://drive.google.com/drive/folders/1FVd05NeUf_zJXD-OlSyinmOi0lOcLdrr?usp=sharing) required for the workshop.

## Breakout Sessions

### Breakout session day 1

Everything we do in AI is driven by the data we have available. Often, this data, whether generated in-house or obtained from literature, is inconsistent and includes many sources of noise. It is important to be aware of the causes for noisy data and how to minimize this before feeding the data into data science tools. It is also helpful to understand the chemical space that the model has trained on and how this relates to the chemical space we are interested in. This understanding will contribute to our confidence in the model predictions before using them for prospective screening.

This activity is split into three parts:

1. A short ice-breaker to get to know your assigned breakout group members.&#x20;
2. A data cleaning activity from a real-world public dataset.
3. A chemical space analysis activity to compare the similarities and differences between chemical libraries.

#### Activity 1: Ice-breaker

Take a few minutes to get to know the members of your breakout group. Take turns to cover the following points about yourself:

1. Your name, research institution and field of research.
2. How you currently use computational tools for your research, if any, and specifically AI-based tools.
3. How you think data science tools can contribute to your research going forward.

Lastly, select a scribe for your breakout group as well as someone else who will provide some feedback on your group’s discussion for the next two tasks during the feedback session.

#### Activity 2: Data Cleaning

<mark style="color:blue;">Task 2.1: Data consistency discussion</mark>

Do some group discussion around the following properties. Why should data have each of them before we train models from it? Think about what problem(s) might result from a dataset that does not have each of these properties.

* Completeness
* Consistency
* Accuracy
* Relevancy

<mark style="color:blue;">Task 2.2: Data cleaning hands-on example.</mark>

During the skills development, we spoke about the need for clean data and some examples of common problems in chemical datasets. Ideally, we need a set of compound structures and corresponding assay outcomes from the same experimental conditions. Now we will identify some examples of data inconsistencies in a dataset from literature.

The Community for Open-Antimicrobial Drug Discovery (CO-ADD) has curated a database of compounds that have been tested for activity against a set of infectious bacteria known as the ESKAPE pathogens. Let’s say we want to create a clean dataset to predict the activity of compounds against the _A.baumannii_ bacteria. Download the dose-response dataset from this [link](https://db.co-add.org/downloads/) and answer the following questions:

1. What is the output of the assay measuring?
2. Find examples of data inconsistencies that cause the dataset to not be complete, consistent, accurate, or relevant for modelling _A.baumannii_ activity. How did you find this issue in the dataset?
3. How could you address each of the points you found in (b) to make the dataset more clean?
4. Were there any other checks you performed where the dataset did not have an issue?

#### Activity 3: Chemical Space Exploration

<mark style="color:blue;">Task 3.1: Chemical space discussion</mark>

In group, discuss and answer the following questions:

1. What do you understand by the term 'chemical space’?
2. How does this differ to the concept of 'drug-like molecules’?
3. Why should the compounds we use to train models and compounds we want predictions for be similar?
4. How might the requirements for our training data vary between the following scenarios?
   1. Virtual screening for novel chemical hits?
   2. Ranking closely related analogues within one chemical series?

<mark style="color:blue;">Task 3.2: Chemical space visualization</mark>

Go to ChEMBL and download the St Jude 3D7 screening set for malaria (ID: CHEMBL730079). Follow the step-by-step instructions in the breakout session slides. Then upload this dataset to the chemical space visualization app (link) and select the ‘MMV Malaria Box’ and ‘Open Source Malaria’ checkboxes from the list of example libraries. Imagine we use the St Jude 3D7 dataset to train a model to screen for new chemical hits in the MMV Malaria Box and Open Source Malaria datasets. Answer the following questions:

1. Do you think that a model that has been trained on the St Jude 3D7 dataset would be more predictive for the MMV Malaria Box or the Open Source Malaria libraries? Why?
2. What steps could you take to make our model more applicable to the library that would be more difficult to make predictions for?
3. How could we improve our model once we’ve experimentally tested the first set of compounds that were selected from the virtual screening?

### Breakout session day 2

A virtual screening cascade allows us to mimic in the computer some of the experimental steps we must do to identify new drug leads. By filtering out molecules with predicted low activities, or undesired side effects, we can lower the cost and time to find new drug candidates. Ideally, we can build virtual screening cascades based off our own data, but for many assays we do not have readily available experimental data. In these situations, we can leverage models developed by third parties and apply them to our problem.

{% hint style="warning" %}
Virtual screening cascades are not meant to substitute experimental testing, but act as a decision-making support tool
{% endhint %}

#### _A.baumannii_ activity prediction <a href="#mmv-malaria-box" id="mmv-malaria-box"></a>

In this activity, we will replicate the work described in Liu et al, 2023, where they build an ML model to identify novel _A.baumannii_ inhibitors and use it to filter the [Drug Repurposing Hub](https://www.broadinstitute.org/drug-repurposing-hub).&#x20;

During the skills development session, we have done a deep dive into ML model building using the _A.baumannii_ model as an example. Now, you have to download the list of compounds available in the Drug Repurposing Hub and continue the "virtual screening" similar to the original author's work. To that end, we suggest running predictions against _A.baumannii_ activity and a few accessory models available through Ersilia to select the best candidates. In short, the steps to follow are:

1. Download the Drug Repurposing Hub data for your group from this [link](https://drive.google.com/drive/folders/1FVd05NeUf_zJXD-OlSyinmOi0lOcLdrr?usp=drive_link).
2. Look at the [Ersilia Model Hub](https://ersilia.io/model-hub) models available online (select Online in the left menu of the website).
3. Select which models would be relevant to your exercise and why, and take notes on model interpretation, expected results, priority level etc.
4. Run predictions for the Drug Repurposing Hub molecules using the selected models.&#x20;
5. Select the best molecule candidates based on your defined filters of activity, ADME properties and other considerations.

{% hint style="info" %}
To simplify the exercise, we have prepared 5 subsets of data from the Drug Repurposing Hub. Each group should use its assigned dataset only
{% endhint %}

#### Models <a href="#models" id="models"></a>

In order to limit the exercise, please limit your screening to the following models:

* A.baumannii Activity: eos3804
* General Antibiotic Activity: eos4e40
* Cardiotoxicity: eos43at
* Synthetic Accessibility: eos9ei3
* ADME properties: eos7d58
* Natural Product Likeness: eos9yui

It is best to select only a few models but really understand how to use them rather than running predictions for all the models but not knowing how to interpret the outcomes.

For each model, think about the following questions:

* What type of model is it (classification or regression)?
* What is the training dataset? (refer to the original publication listed above)
* What is the interpretation of the model outcome?
* What cut-off, if any, we should use for that particular model?

In addition, think about the following concepts:

* Does the outcome of the model make sense? If it does not make sense, perhaps we have the wrong interpretation of the model output.
* Is the cut-off I have selected too stringent (i.e, I am losing too many molecules and I should be more permissive?)
* Is this model very relevant for the current dataset (i.e, is malaria activity equally important as natural product likeness?)

**Molecule selection**

Use the predicted values to select the 10 molecules that you would take for experimental testing if you had to choose. To that end, you can think of:

* What are the most important activities you want to optimize
* What are strict no-go points
* What are activities that are easiest to optimize at lead stage

Finally, prepare a short presentation for the rest of the participants. This should cover:

* Which models did you choose and why
* What selection strategy did you decide
* Which were your selected molecules

{% hint style="success" %}
Extension: if you finish the proposed activity, have a look at what else is available in the Ersilia Model Hub and what would you like to see deployed online!
{% endhint %}

### Breakout session day 3

Sampling the chemical space is a fundamental task in computational drug discovery. Recently, generative AI methods have significantly expanded our sampling capabilities. As a result, traditionally laborious steps, such as hit-to-lead optimization, can now be supported by computational tools that generate new chemical matter as starting points for further analysis. In this session, we will experiment with some of these tools to sample the chemical space around a seed molecule.

#### Step 1: Sample the chemical space around a seed molecule

You will be presented with four hits obtained from an experimental screening against _Burkholderia cenocepacia_ (growth inhibition assay). These will serve as your starting points. As a team, **agree on one hit (seed) molecule for further analysis**. There is no right or wrong answer.

Next, you will have the opportunity to experiment with different “samplers” from the Ersilia Model Hub. Read the model descriptions and discuss the following:

* Some samplers perform a similarity search against a large chemical library. Can you identify these samplers? Do you expect them to give deterministic results?
* Other samplers generate genuinely new molecules. Of these, two are fragment-based and two are deep-learning (AI) based. Can you classify them?
* What are the expected advantages and disadvantages of each type of model?
* Some samplers may not yield results. What might be the reason for this?

#### Step 2: Create a wishlist of molecules

The best way to understand the behavior of sampling methods is to run them and explore the results, both visually and quantitatively. Feel free to **sample compounds multiple times, trying different methods**. There are filtering options in the app, as well as a simple display of the molecular structures.

Ultimately, you should **create a wishlist of 100-500 molecules**, ideally obtained using multiple models. To create this wishlist, consider the following questions:

* What is a reasonable range of Tanimoto similarities to the seed molecule?
* Of the auxiliary properties (e.g., molecular weight, logP, and QED), which is the most relevant at this stage?

{% hint style="info" %}
Pro tip: open an Excel spreadsheet and copy-paste your molecules of interest them. Once you are satifsfied with them, copy-paste them in the relevant box in the app.
{% endhint %}

#### Step 3: Predict properties of your molecule wishlist

The Ersilia Model Hub contains a variety of activity and property prediction models that can be used to further assess your wishlist. In this case, we have selected **three activity prediction models** ([eos4e40](https://github.com/ersilia-os/eos4e40), [eos5xng](https://github.com/ersilia-os/eos5xng), and [eos9f6t](https://github.com/ersilia-os/eos9f6t)) and a multi-output **ADMET model** ([eos7d58](https://github.com/ersilia-os/eos7d58)). Discuss the following with your team:

* Which activity prediction model is best suited for this task?
* Which 5 ADMET properties would you prioritize? There is no right or wrong answer.
* What is more important at this stage: activity or ADMET properties?
* ADMET values will be given in percentiles. How do you interpret these values?

Once you have made your decision, you can calculate your desired enpoints for all the molecules in your wishlist.

Explore results on two levels:

1. **General level:** Assess whether wishlist molecules have high predicted activities (as desirable) and acceptable ADMET properties.
2. **Individual level:** Select 3-5 molecules and discuss whether they are good leads.

{% hint style="success" %}
Send your selected 3-5 compounds to [miquel@ersilia.io](mailto:miquel@ersilia.io). The title of your email should be "AI2050 Group N". Just copy-paste the SMILES strings of the molecules in the e-mail text.
{% endhint %}

#### Final step: Prepare a presentation

Prepare a short presentation for the rest of the participants. The presentation should cover the following:

* What seed molecule did you choose and why?
* Which sampler models did you find most useful?
* What was your rationale for assembling your wishlist?
* Which lead molecules did you select and why?

