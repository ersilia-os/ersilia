---
description: >-
  This page describes the contribution guidelines for the interns interested in
  participating in the Outreachy round from December 2022 to March 202
---

# Outreachy Winter 2023

## Internship project

The Ersilia Model Hub is an open source platform of ready-to-use AI/ML models for biomedical research. With it, scientists can browse a collection of models, select the ones relevant to their research and run predictions! For example, to predict whether a molecule will be active against Malaria.

The internship project will be focused on increasing the collection of available models in the Hub.&#x20;

Specially important to Ersilia is the contribution period. We list below a number of tasks that must be completed in order. <mark style="color:red;">We will not judge interns on how many tasks they can complete, but on the quality of each contribution and the interest to learn, participate and help others in the community.</mark> The mentors are there and willing to help, but they also have a dedicated timeslot for reviewing contributions, please be patient if it takes a few hours to get back to you.

{% hint style="warning" %}
Please note that this internship requires knowledge of Python programming language
{% endhint %}

## Signing up for Outreachy

Interested applicants must be accepted in the Outreachy internship program. Please go to Outreachy's [website](https://outreachy.org) and be sure to fill in the application according to the timelines.

{% hint style="danger" %}
Ersilia can only accept interns that have been approved by Outreachy and that comply with the necessary requirements. Please check your availability for the internship period fulfills Outreachy's requirements.
{% endhint %}

## Contribution period

The contribution period runs from <mark style="color:blue;">October 6th</mark> to <mark style="color:blue;">November 4th</mark>. During this time, interested applicants are welcome to contribute to Ersilia's project following the guidelines in this document.

### 1. Get to know the community

We will be having tons of interactions during the contribution period, so the best is to get to know your peers, your mentors and a bit more about Ersilia.

#### Slack

We use Slack as our main communication platform, both for the contribution period and afterwards to work with the selected interns. If you have never used this tool, don't worry, is quite intuitive!

1. [Sign up](https://join.slack.com/t/ersilia-outreachy/shared\_invite/zt-1h2u8vbxh-fkyxlOJE20Bx538hQvpC\~A) using your preferred email and name.
2. Introduce yourself in the #general channel. The general channel is for random questions, interactions with other fellow contributors, writing tips and suggestions...
3. Use the dedicated channels for questions about specific topics. For example, the #colab channel will be used to discuss about Google Colaboratory issues.&#x20;
4. Contribute to your peers questions, this is about helping each other and we really value interns who work with the community.

{% hint style="info" %}
We try to work as openly as possible, we encourage all contributors to post in the open channels rather than private conversations.
{% endhint %}

{% hint style="success" %}
Please use a Slack name that is easy to identify with your GitHub handle to make it easy for mentors to review contributions and tag people.
{% endhint %}

#### GitHub

We will be working on the Ersilia Model Hub main repository, which is hosted on [GitHub](https://github.com/ersilia-os/ersilia). You can start by:

* 📖 Getting familiar with the repository structure
* 🐛 Checking the issues to see what has the community been working on
* 👀 Watching the repository to receive notifications if you are mentioned&#x20;
* ⭐ Starring the repository if you like the work Ersilia is doing!

{% hint style="danger" %}
Ersilia is an Open Source community with active contributing members. Please respect the work of others and specially if an issue is assigned to someone else, give them the time and space to work on it.
{% endhint %}

{% hint style="info" %}
We will be using GitHub issues a lot, so if you have never worked with GitHub before, make sure you understand how issues work!
{% endhint %}

#### Community Call

We will host an open community call on <mark style="color:blue;">Wednesday 12th of October</mark> at 17:00 SAST. During the call we will introduce a bit more about Ersilia's work, go over the contribution period tasks and answer any questions you might have!

The link to the call will be shared via Slack.

We apologize in advance if the call is in a difficult time-zone for you, attendance is not compulsory and will not be taken into account when selecting the interns.&#x20;

#### Code of Conduct

Ersilia is adhered to the [Contributor Covenant Code of Conduct](../../about-us/code-of-conduct.md). Any breaches of the code of conduct, specially harassment or lack of respect for fellow contributors, will have severe consequences.

### 2. Install the Ersilia Model Hub

We will be using the Ersilia Model Hub throughout the internship. The software is in beta testing and therefore you might encounter errors while running it, don't worry this is why we are here!

Please follow the installation [instructions](../../ersilia-model-hub/installation.md). If you have a UNIX machine (Linux or MacOS) you can install Ersilia directly. If you are using a windows machine you will need a Virtual Machine or a Windows Subsystem Linux (WSL).

{% hint style="info" %}
For Windows users, we recommend using a WSL with Visual Studio Code to access it.
{% endhint %}

{% hint style="danger" %}
A common mistake is to forget the installation of Git-LFS, which is required for many models. Please do so!
{% endhint %}

#### Testing that Ersilia works

We will first make sure ersilia works by running the following commands:

```
ersilia --help #this should output the command options for ersilia
```

Once we are sure ersilia is recognised in the CLI, we will test a very simple model

```
ersilia -v fetch eos3b5e
ersilia serve eos3b5e
ersilia -v api calculate -i "CCCC"
```

This is calculating the molecular weight of the molecules, the output should be printed in your CLI and look like:

```
{
    "input": {
        "key": "IJDNQMDRQITEOD-UHFFFAOYSA-N",
        "input": "CCCC",
        "text": "CCCC"
    },
    "output": {
        "mw": 58.123999999999995
    }
}
```

{% hint style="danger" %}
These tests do not work, what now?! Open an issue on GitHub, indicating on the title Ersilia installation problem and giving a full description of the errors you get. The mentors will answer as soon as possible.
{% endhint %}

{% hint style="success" %}
These work just fine! Perfect, tag the mentor Gemma on the #stage1-contributors Slack channel so she can assign you models for the next step (see below)
{% endhint %}

### 3. Run models on the Ersilia CLI

Your second task for the contribution period is to help us debug any issues when running models with the Ersilia Command Line Interface. You can find more information about it on the [Model Usage](../../ersilia-model-hub/antibiotic-activity-prediction/) guide.

{% hint style="warning" %}
Use the -v (verbose) in all commands to see the error outputs.

`ersilia -v fetch ...`

`ersilia -v serve ...`

`ersilia -v api ...`
{% endhint %}

To keep track of the models we have tested, please add your information to the [shared excel file](https://docs.google.com/spreadsheets/d/1IsWv63Krt3g6QcnTTK9OFA-5IBKVp2y7HD5DgCYLoGc/edit?usp=sharing). In short:

1. Use one of the empty \<username> column fields
2. Add the system you are using on the first row.
3. Change the cells of your column corresponding to the models you have tested in:
   1. <mark style="color:green;background-color:green;">Green</mark><mark style="color:green;">,</mark> if they run without a problem
   2. <mark style="color:red;background-color:red;">Red</mark>, if there was an issue
4. If you run the model <mark style="color:green;">successfully:</mark>&#x20;
   1. Write fetch and predict times on the excel cell.
   2. Go on to the next model, until you have tested the models assigned to you.
5. If you encountered a problem:
   1. Go to Ersilia's GitHub repository and check if there is a Bug already for this model. If there is, add your error on the same thread.
   2. If not, open a Bug issue (use the provided template) and add the model identifier (eosxxxx) as issue title.
   3.  Add the log of the error, which can be obtained by adding the following code snippet at the end of the command&#x20;

       ```
       ersilia -v fetch modelname > my.log 2>&1
       ```
   4. Let's work on debugging this issue together before moving on!

We have provided an already prepared list of [test molecules](https://drive.google.com/file/d/14bZ1mcCv0jD-Wc-\_qwz-uxMGrJaenaup/view?usp=sharing) in .csv format. Download it and use it for model testing. Store the output of the model in a .csv file that must contain in the name the <mark style="color:orange;">model identifier</mark> and the <mark style="color:orange;">date</mark> when the prediction was done. Upload the output .csv file in the [shared folder.](https://drive.google.com/drive/folders/1h43ndbpPWkz2rIHHrb0H-5rgG8cpu803?usp=sharing)

{% hint style="warning" %}
Do not overwrite your fellow's additions to the excel file
{% endhint %}

{% hint style="danger" %}
To help the mentors keep track of your progress, add all relevant information to the GitHub issue
{% endhint %}

### 4. Run models on Google Colab

Once you have successfully tested 5 models, go on and try to repeat the exercise in Google Colaboratory. You can also test one model in the CLI and in Google Colab in parallel. You can use the[ Colaboratory Template Notebook](https://github.com/ersilia-os/ersilia/blob/master/notebooks/ErsiliaOnColaboratory.ipynb) from Ersilia's repository.

[Google Colaboratory](../../training-materials/google-colaboratory.md) uses Google's servers (Linux machines) to run the code, which is very convenient to bypass installation issues. To run Google Colab, you only need a Google Account. If you do not have and do not wish to open a Google Account, we can skip this step.All

{% hint style="info" %}
All the model requirements stay within Google's free tier
{% endhint %}

### 5. Train a new model and incorporate it in the Hub

Once you have:

* Successfully installed the Ersilia Model Hub in your computer
* Successfully run predictions for at least, 3 models using the command line AND debugged any issues you might have.
* Successfully run the same five models in Google Colaboratory using the Ersilia PyPi package

We are ready to continue onto the next stage of the contribution period 🎉

For this period, we will use Ersilia's automated ML modelling packages:

* [lazy-qsar](https://github.com/ersilia-os/lazy-qsar): a library to build rapid QSAR models
* [zairachem](https://github.com/ersilia-os/zaira-chem): an end-to-end QSAR modelling library with excellent performance rates

We will leverage the datasets from the excellent initiative Therapeutics Data Commons ([TDC](https://tdcommons.ai/)). You can read more about it in its associated [publication](https://www.nature.com/articles/s41589-022-01131-2). TDC has prepared biomedical-related datasets for ML modelling, and provides benchmarks of performance. We will use those to test our automated ML libraries and add the resulting models in the Ersilia Model Hub.

For this part, please do not launch directly onto modelling, wait for the mentor's revision of your previous work:

1. Once you have completed all the above steps, mention the mentors on the Slack channel #stage1-contributions and explain all the steps done on the first phase. Add your GitHub handle in the issue!
2. The mentors will open a GitHub issue with the specific modelling exercise and assign it to you.
3. All the conversations related to developing the models for the dataset assigned to you will be discussed in the GitHub issue

Once the model is completed, we will add it to the Ersilia Model Hub. Read more about it in the [guidelines](../incorporate-models.md).
