---
description: >-
  This page describes the contribution guidelines for the interns interested in
  participating in the Outreachy round from May 2023 to August 2023
---

# Outreachy Summer 2023

## Internship project

The Ersilia Model Hub is an open source platform of ready-to-use AI/ML models for biomedical research. With it, scientists can browse a collection of models, select the ones relevant to their research and run predictions! For example, to predict whether a molecule will be active against Malaria. The internship project will be focused on increasing the collection of available models in the Hub.&#x20;

Specially important to Ersilia is the contribution period. We list below a number of tasks that must be completed in order. <mark style="color:red;">We will not judge interns on how many tasks they can complete, but on the quality of each contribution and the interest to learn, participate and help others in the community.</mark> The mentors are there and willing to help, but they also have a dedicated timeslot for reviewing contributions, please be patient if it takes a few hours to get back to you. **Almost all the information you will need to successfully complete the contribution period is enclosed in this document, please read it all before asking questions.**

{% hint style="warning" %}
Please note that this internship requires knowledge of Python programming language.

Also, all text written by interns will be revised to ensure no plagiarism and no use of AI writing tools has been involved. While AI tools can be helpful in certain circumstances, here we need to learn about you in your own words. Using Chat-GPT or others will immediately disqualify you from further participating in the program.
{% endhint %}

## Signing up for Outreachy

Interested applicants must be accepted in the Outreachy internship program. Please go to Outreachy's [website](https://outreachy.org) and be sure to fill in the application according to the timelines.

{% hint style="danger" %}
Ersilia can only accept interns that have been approved by Outreachy and that comply with the necessary requirements. Please check your availability for the internship period fulfills Outreachy's requirements.
{% endhint %}

## Contribution period

The contribution period runs from <mark style="color:blue;">March 6th</mark> to <mark style="color:blue;">March 31st</mark>. During this time, interested applicants are welcome to contribute to Ersilia's project following the guidelines in this document.

The contribution period is organised in 4 weeks. Each week has a set of specific goals defined, with the objective that mentors can evaluate the intern's experience, interest in the community and team-building work. Once the week's objectives have been met, please focus on:&#x20;

* Improving your contribution (there is always more publications to read, better bug reports to be written etc)&#x20;
* Helping out other contributors (we really value group work)

{% hint style="info" %}
Link this issue to your outreachy application
{% endhint %}

### 📆 WEEK 1: Get to know the community

The first week is focused on getting to know the Ersilia community, our mission and how we work to achieve it. We will be having tons of interactions during the contribution period, so the best is to get to know your peers, your mentors and a bit more about Ersilia.

#### Task 1: Join in the communication channels

**Slack:** we use Slack as our main communication platform, both for the contribution period and afterwards to work with the selected interns. If you have never used this tool, don't worry, is quite intuitive!

1. [Sign up](https://join.slack.com/t/ersilia-outreachy/shared\_invite/zt-1h2u8vbxh-fkyxlOJE20Bx538hQvpC\~A) using your preferred email and name.
2. Introduce yourself in the #general channel. The general channel is for random questions, interactions with other fellow contributors, writing tips and suggestions...
3. Use the dedicated channels for questions about specific topics.&#x20;
4. Contribute to your peers questions, this is about helping each other and we really value interns who work with the community.

{% hint style="info" %}
We try to work as openly as possible, we encourage all contributors to post in the open channels rather than private conversations.
{% endhint %}

{% hint style="success" %}
Please use a Slack name that is easy to identify with your GitHub handle to make it easy for mentors to review contributions and tag people.
{% endhint %}

**GitHub**: we will be working on the Ersilia Model Hub main repository, which is hosted on [GitHub](https://github.com/ersilia-os/ersilia). You can start by:

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

**Community Call:** we will host an open community call (TBD). During the call we will introduce a bit more about Ersilia's work, go over the contribution period tasks and answer any questions you might have! The link to the call will be shared via Slack. We apologize in advance if the call is in a difficult time-zone for you, attendance is not compulsory and will not be taken into account when selecting the interns.&#x20;

**Code of Conduct:** Ersilia is adhered to the [Contributor Covenant Code of Conduct](../../about-us/code-of-conduct.md). Any breaches of the code of conduct, specially harassment or lack of respect for fellow contributors, will have severe consequences.

#### Task 2: Open an issue

Each intern will have to open an issue on the Ersilia Model Hub repository using the Outreachy Summer 23 [issue template](https://github.com/ersilia-os/ersilia/issues/new/choose). Use this issue as your personal thread of the contribution period. Everything you want the mentors to review, comment or take into account for the contribution period should be in that issue. You can also tag people to ask for help, debug together...

The issue has a set of tasks, mark them as complete as you go!

#### Task 3: Install the Ersilia Model Hub

We will be using the Ersilia Model Hub throughout the internship. The software is in beta testing and therefore you might encounter errors while running it, don't worry this is why we are here!

Please follow the installation [instructions](../../ersilia-model-hub/installation.md). If you have a UNIX machine (Linux or MacOS) you can install Ersilia directly. If you are using a windows machine you will need a Virtual Machine or a Windows Subsystem Linux (WSL).

{% hint style="info" %}
For Windows users, we recommend using a WSL with Visual Studio Code to access it.
{% endhint %}

{% hint style="danger" %}
A common mistake is to forget the installation of Git-LFS, which is required for many models. Please do so!
{% endhint %}

1. Testing that Ersilia works: we will first make sure ersilia works by running the following commands:

```
ersilia --help #this should output the command options for ersilia
```

2. Once we are sure ersilia is recognised in the CLI, we will test a very simple model

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
These tests do not work, what now?! Write down the challenges you are facing in your github issue, and ask for support to your peers through the Slack channel.
{% endhint %}

{% hint style="success" %}
These work just fine! Perfect, report the tasks as completed in your specific GitHub issue and help your peers achieve it as well!
{% endhint %}

#### Task 4: Motivation statement

Write, in a thread in your issue, your motivation for joinign Outreachy and, in particular, why are you interested in working at Ersilia. A good motivation letter will explain your current skills, your reasons to work in Ersilia's project, how this would advance your career and what are your plans during and after the internship.

#### Task 5: open an application to Ersilia

If what you have seen and learnt in this first week is appealing and you want to continue working with us, open an application to Ersilia through the Outreachy website, and link the issue you have opened to track your contributions

### 📆 WEEK 2: Run an ML model

The task of this week is to successfully install an ML model in your system (using third-party code). The goal of this week is to:

* Demonstrate your python knowledge
* Practise following documentation and implementing third party code
* Working with dependencies and conda environments

This task is very similar to what you will work on during the internship; finding open source code, testing it, debugging if necessary and finally, implementing it within Ersilia.

#### Task 1: select a model

We propose one of the following models to be implemented. Please select one and explain in the GitHub issue why you have selected it (interest in the application, in the ML algorithm used...?)

* NCATS Rat Liver Microsomal Stability: [https://github.com/ncats/ncats-adme](https://github.com/ncats/ncats-adme)
* Plasma Protein Binding (IDL-PPBopt): [https://github.com/Louchaofeng/IDL-PPBopt](https://github.com/Louchaofeng/IDL-PPBopt)
* SARS-CoV2 activity (Image Mol): [https://github.com/HongxinXiang/ImageMol](https://github.com/HongxinXiang/ImageMol)
* STOUT (SMILES to IUPAC):: [https://github.com/Kohulan/Smiles-TO-iUpac-Translator](https://github.com/Kohulan/Smiles-TO-iUpac-Translator)

{% hint style="warning" %}
These are Open Source models developed by academics mostly, and they have nothing to do with Ersilia or the Outreachy program, so please refrain from contacting the original authors without asking one of the Ersilia mentors first.
{% endhint %}

#### Task 2: install the model

Follow the instructions on the specific model repository to install and run the model in your system. Report the issues you find along the process and how you solve them in the GitHub issue.

A good task here will detail each and every one of the steps taken to successfully run and implement the model in your local machine. If you simply mark the task as complete, the mentors will not be able to evaluate your knowledge.

#### Task 3: run predictions

Run predictions for the [Essential Medicines List](https://raw.githubusercontent.com/ersilia-os/ersilia/master/notebooks/eml\_canonical.csv) provided in our repository, post the result of these predictions (using the original source code) in the GitHub thread and explain what they mean.

{% hint style="info" %}
Showing that you are not only able to run code but actually understand the outputs of the ML models is critical for the intership
{% endhint %}

#### Task 4: compare results with the Ersilia models

The models we suggest are models that are already incorporated in the [Ersilia Model Hub](https://ersilia.io/model-hub). Look for them, fetch them, run predictions and compare the output of the original model with the output of the original code.

### 📆 WEEK 3: Literature search

A big part of what we do at Ersilia is to screen the scientific literature in search of new models and datasets of interest to our community. We are always looking for models that can help speed up drug discovery against infectious and neglected diseases.

In this week, focus on diving into the scientific literature, trying to find studies of interest to our community. To that end, we suggest using websites like PubMed, Google Scholar, BioRxiv and Papers with Code.&#x20;

To get inspired, check the [Contribute models](../incorporate-models.md) section, look at [existing models](https://ersilia-io/model-hub) in the Hub and our list of [pending models](https://airtable.com/shrTpe45mLKqaHXsc) to incorporate.

#### Task 1: a first model suggestion

Find one publication that describes an open source ML model that could be of interest to Ersilia (activity against a specific pathogen, cytotoxicity, side effects...) and link it in the thread.

Explain what the model does (in your own words), why it would be relevant to Ersilia and how would you implement it (look at its code and whether it is ready to be used, e.g are the model checkpoints provided, is the underlying data available?)

{% hint style="danger" %}
Do not open a nw model request issue, or add it to our backlog. Only do so if mentors review the model suggestion and ask you to do it.
{% endhint %}

#### Tasks 2 and 3:

Continue your model search and suggest up to three studies/models/publications of interest, each one with detailed explanations.

### 📆 WEEK 4: submit your final application

Focus the last week ONLY in writing your final application to Outreachy. Mentors will not revise any contribution for the last week, only final applications.

We provide this template of an 8-week internship at Ersilia. Use it and adapt it to your internship plans (reminder: Outreachy lasts 12 weeks).

{% hint style="warning" %}
Final applications must be submitted through the Outreachy website on time. We will not be able to provide help or support for last minute internet connection problems, late submissions and other issues. Please make sure you fill it in with time.
{% endhint %}
