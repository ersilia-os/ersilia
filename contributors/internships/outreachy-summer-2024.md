---
description: >-
  Please find here the guidelines for the Outreachy contribution period running
  from 4th March to 2nd April 2024.
---

# Outreachy Summer 2024

## Internship Project

The Ersilia Model Hub is an open source platform of ready-to-use AI/ML models for biomedical research. With it, scientists can browse a collection of models, select the ones relevant to their research and run predictions! For example, to predict whether a molecule will be active against Malaria. The internship project will be focused on increasing the collection of available models in the Hub.

Specially important to Ersilia is the contribution period. We list below a number of tasks that must be completed in order. <mark style="color:red;">We will not judge interns on how many tasks they can complete, but on the quality of each contribution and the interest to learn, participate and help others in the community.</mark> The mentors are there and willing to help, but they also have a dedicated time slot for reviewing contributions, please be patient if it takes a few hours to get back to you. **Almost all the information you will need to successfully complete the contribution period is enclosed in this document, please read it all before asking questions.**

{% hint style="warning" %}
Please note that this internship requires knowledge of Python programming language, and a basic understanding of fundamental data science concepts.
{% endhint %}

{% hint style="danger" %}
Please note that we have a zero tolerance policy for plagiarism. All text written by interns will be revised to ensure no plagiarism and no use of AI writing tools has been involved. While AI tools can be helpful in certain circumstances, here we need to learn about you in your own words. Using Chat-GPT or others for writing letters of interest, summaries of publications or similar will immediately disqualify you from further participating in the program.
{% endhint %}

## Signing up for Outreachy

Interested applicants must be accepted in the Outreachy internship program. Please go to Outreachy's [website](https://outreachy.org) and be sure to fill in the application according to timelines.

{% hint style="danger" %}
Ersilia can only accept interns that have been approved by Outreachy and that comply with the necessary requirements. Please check your availability for the internship period fulfills Outreachy's requirements.
{% endhint %}

## Contribution Period

The contribution period runs from <mark style="color:blue;">March 4th</mark> to <mark style="color:blue;">April 2nd</mark>. During this time, interested applicants are welcome to contribute to Ersilia's project following the guidelines in this document.

The contribution period is organised in 4 weeks. Each week has a set of specific goals defined, with the objective that mentors can evaluate the intern's experience, interest in the community and team-building work. Once the week's objectives have been met, please focus on:

* Improving your contribution (there is always more publications to read, better bug reports to be written etc)
* Helping out other contributors (we really value group work)

We will be using GitHub issues to track the work of each contributor. Please open a new issue on our Ersilia [repository](https://github.com/ersilia-os/ersilia), choosing the Outreachy Summer 2024 [template](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=internship\&projects=\&template=outreachy_s24.yml\&title=%E2%9C%8D%EF%B8%8F+Contribution+period%3A+%3Cyour_name%3E+). You can check the tasks while you complete them.

{% hint style="info" %}
Link this issue to your Outreachy application to submit the first contribution.
{% endhint %}

### 📆 WEEK 1: Get to know the community

The first week is focused on getting to know the Ersilia community, our mission and how we work to achieve it. We will be having tons of interactions during the contribution period, so the best is to get to know your peers, your mentors and a bit more about Ersilia.

#### Task 1: Join in the communication channels

**Slack:** we use Slack as our main communication platform, both for the contribution period and afterwards to work with the selected interns. If you have never used this tool, don't worry, is quite intuitive!

1. [Sign up](https://join.slack.com/t/ersilia-outreachy-s24/shared_invite/zt-2e1h1cjs5-lGyrZPprqlpL~qEWgeoCWw) using your preferred email and name.
2. Introduce yourself in the #general channel. The general channel is for introductions, announcements, random questions, interactions with other fellow contributors, writing tips and suggestions.
3. Use the dedicated channels for questions about specific topics.
4. Contribute to your peers' questions, this is about helping each other and we really value interns who work with the community.

{% hint style="info" %}
We try to work as openly as possible, we encourage all contributors to post in the open channels rather than private conversations.
{% endhint %}

{% hint style="info" %}
Please use a Slack name that is easy to identify with your GitHub handle to make it easy for mentors to review contributions and tag people.
{% endhint %}

**GitHub**: we will be working on the Ersilia Model Hub main repository, which is hosted on [GitHub](https://github.com/ersilia-os/ersilia). You can start by:

* 📖 Getting familiar with the repository structure
* 🐛 Checking the issues to see what has the community been working on
* 👀 Watching the repository to receive notifications if you are mentioned
* ⭐ Starring the repository if you like the work Ersilia is doing!

{% hint style="danger" %}
Ersilia is an Open Source community with active contributing members. Please respect the work of others and specially if an issue is assigned to someone else, give them the time and space to work on it.
{% endhint %}

{% hint style="info" %}
We will be using GitHub issues a lot, so if you have never worked with GitHub before, make sure you understand how issues work!
{% endhint %}

**Community call:** we will hold a community call on **March 8th 10:00 am CET** to go over the contribution period tasks and answer any questions you might have! The link will be shared via Slack. Attendance is not compulsory, we have tried to find a time that is acceptable for most time zones, we apologize in advance if it means an early start or late end of your day.

**Code of Conduct:** Ersilia is adhered to the [Contributor Covenant Code of Conduct](../../about-us/code-of-conduct.md). Any breaches of the code of conduct, specially harassment or lack of respect for fellow contributors, will mean disqualification as an applicant.

#### Task 2: Open an issue

Each intern will have to open an issue on the Ersilia Model Hub repository using the Outreachy Summer 24 issue [template](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=internship\&projects=\&template=outreachy_s24.yml\&title=%E2%9C%8D%EF%B8%8F+Contribution+period%3A+%3Cyour_name%3E+). Use this issue as your personal thread of the contribution period. Everything you want the mentors to review, comment or take into account for the contribution period should be in that issue. You can also tag people to ask for help, debug together...

The issue has a set of tasks, mark them as complete as you go!

#### Task 3: Install the Ersilia Model Hub

We will be using the Ersilia Model Hub throughout the internship. The software is in beta testing and therefore you might encounter errors while running it, don't worry this is why we are here!

Please follow the installation [instructions](broken-reference). If you have a UNIX machine (Linux or MacOS) you can install Ersilia directly. If you are using a windows machine you will need a Virtual Machine or a Windows Subsystem Linux (WSL).

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
ersilia -v api run -i "CCCC"
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
These tests do not work, what now?! Write down the challenges you are facing in your GitHub issue, and ask for support to your peers through the Slack channel.
{% endhint %}

{% hint style="success" %}
These work just fine! Perfect, report the tasks as completed in your specific GitHub issue and help your peers achieve it as well!
{% endhint %}

#### Task 4: Testing Ersilia with Docker

Ersilia models are supported to be fetched from a whole host of places, namely, S3 buckets, GitHub, DockerHub, and even from a local repository of the model!

To ensure model dependencies are self contained, Ersilia models are "dockerized", and Ersilia fetches them through DockerHub by default, if you have Docker installed. To complete this task, make sure you have Docker installed or install it from [here.](https://www.docker.com/get-started/)

1. Pull a model image. Here we use another simple model from the hub

```
docker pull ersiliaos/eos4wt0:latest
```

2. Activate the environment where you have installed ersilia, and test this model that you have just fetched from DockerHub:

```
ersilia serve eos4wt0 # Notice that you don't have to fetch it through ersilia here.
ersilia -v api run -i "CCCC"
```

This generates the Morgan Fingerprints for a molecule and the output should be printed in your shell like this:

```
{
    "input": {
        "key": "IJDNQMDRQITEOD-UHFFFAOYSA-N",
        "input": "CCCC",
        "text": "CCCC"
    },
    "output": {
        "outcome": [
            0.0,
            0.0,
            0.0,
            0.0,
            ...,
        ]
    }
}
```

#### Task 5: Motivation statement

Write, in a thread in your issue, your motivation for joining Outreachy and, in particular, why are you interested in working at Ersilia. A good motivation letter will explain your current skills that are relevant to Ersilia, your reasons to work in Ersilia's project, how this would advance your career and what are your plans during and after the internship.

#### Task 5: open an application to Ersilia

If what you have seen and learnt in this first week is appealing and you want to continue working with us, open an application to Ersilia through the Outreachy website, and link the issue you have opened to track your contributions.

{% hint style="warning" %}
We will not be providing further feedback on your task if you do not record the first contribution on the Outreachy website. This is to ensure we focus our limited capacity for support of those candidates who want to make a final application to the organisation
{% endhint %}

### 📆 WEEK 2: Get Familiar with Machine Learning for Chemistry

The task for this week is to get familiar with Machine Learning for chemistry data, which is going to be the focus of the internship. The Ersilia Model Hub contains about 150 models, which predict several properties of small molecules (bioactivity against pathogens, ADME properties, toxicity…). We want to make sure those models are accurate and reproducible. To that end, we focus on validating the models in three steps:

* **T1** Model bias (i.e: models giving very high values or low values): to check that, we only need to run predictions for a list of 1000 diverse molecules in each model and plot the results in a scatter plot.
* **T2** Reproducibility: can we reproduce the exact values / a figure / that authors obtained when training the model in the first place? This means we need to read the publication and identify for example a compound identified using that model and check that we obtain the same values.
* **T3** Performance: can we check if the model gives accurate results in external datasets? This is more time consuming and will be done by identifying a public dataset that has not been used in model training, and running predictions to build AUROC curves - to simplify reports, we will only focus on AUROC or R2 as metrics now.

For week 2, we will focus on Task 1 and Task 2.&#x20;

Task 1:

1. Pick a model from the list:&#x20;
   1. hERG Models:
      1. [eos2ta5](https://github.com/ersilia-os/eos2ta5)
      2. [eos30f3](https://github.com/ersilia-os/eos30f3)
      3. [eos30gr](https://github.com/ersilia-os/eos30gr)
      4. [eos4tcc](https://github.com/ersilia-os/eos4tcc)
   2. ADME Models:
      1. [eos74bo](https://github.com/ersilia-os/eos74bo)
      2. [eos31ve](https://github.com/ersilia-os/eos31ve)
      3. [eos9tyg](https://github.com/ersilia-os/eos9tyg)
      4. [eos6oli](https://github.com/ersilia-os/eos6oli)
2. Using [this repository](https://github.com/GemmaTuron/model-validation-example) as an example, create a repository in your GitHub account with the appropriate structure (for example, the repository will contain a README file that will explain what is in there, steps to run/reproduce your work, the necessary datasets and folders, licenses and gitignore files)
3. Download and run the selected model from Ersilia and make sure it works
4. Select a list of 1000 molecules from public repositories and make sure they are represented as standard SMILES&#x20;
5. Run predictions for the 1000 molecules, create the necessary plots and explain the results you are obtaining

_Hints:_

1. _Try to understand what the problem that the model is trying to solve, and how this problem is formulated? For example, is the model trying to generate a score, or a probability._
2. _Work on Jupyter notebooks to make it easy for the mentors to review the work._

Once Task 1 is complete, ask for review from the mentors and incorporate any feedback you are given.&#x20;

Task 2:

1. Read the scientific publication and identify a result you could reproduce from the paper. _Hint: explain in your GitHub issue for which dataset are you going to reproduce the results given in the publication._
2. Try to implement the model as described by the authors (not from the Ersilia Model Hub) and reproduce the results selected from the paper. Make sure to include the results evaluation in your repository.
3. Check that the model provides the same results when running via the Ersilia Model Hub.

Once tasks 1 and 2 are complete, move onto tasks for Week 3, and yes you can start even if the week hasn't started yet! :)&#x20;

### 📆 WEEK 3: Validate a Model in the Wild

This is a continuation from the efforts of Week 2. During this week, you will work with the mentors to find a dataset for which:

* We have experimental results on the activity we are trying to predict (for example, malaria IC50)
* There is no "data leakage", ie, molecules in this dataset are not in the training set of the model. This is a very important check!

Once the dataset is curated and cleaned, we can run predictions and create performance metrics like AUROC curves.

### 📆 WEEK 4: Submit your Final Application

Focus the last week ONLY in writing your final application to Outreachy. Mentors will not revise any contribution for the last week, only final applications, to ensure we can provide feedback on them.

We provide [this template](https://docs.google.com/document/d/1-ubLY_LlSApClrjYzmoaR0M-3LTykRi71xBEi4wvO8Q/edit) of a former 8-week internship at Ersilia as an example. Use it and adapt it to your internship plans (reminder: Outreachy lasts 12 weeks).

{% hint style="warning" %}
Final applications must be submitted through the Outreachy website on time. We will not be able to provide help or support for last minute internet connection problems, late submissions and other issues. Please make sure you fill it in with time.
{% endhint %}
