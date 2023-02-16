---
description: >-
  This page describes the contribution guidelines for the interns interested in
  participating in the Outreachy round from May 2023 to August 2023
---

# Outreachy Summer 2023

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

We will host an open community call (TBD). During the call we will introduce a bit more about Ersilia's work, go over the contribution period tasks and answer any questions you might have!

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

### 3. Next steps

This section is under development. More details on specific tasks to be tackled by interns will be provided soon ⚠️
