---
description: >-
  Please find here the guidelines for the Outreachy contribution period running
  from 1st to 29th October
---

# Outreachy Winter 2024

## Internship Project

The Ersilia Model Hub is an open source platform of ready-to-use AI/ML models for biomedical research. With it, scientists can browse a collection of models, select the ones relevant to their research and run predictions! For example, to predict whether a molecule will be active against Malaria. The internship project will be focused on increasing the collection of available models in the Hub.

Specially important to Ersilia is the contribution period. We list below a number of tasks that must be completed in order. <mark style="color:red;">We will not judge interns on how many tasks they can complete, but on the quality of each contribution and the interest to learn, participate and help others in the community.</mark> The mentors are there and willing to help, but they also have a dedicated time slot for reviewing contributions, please be patient if it takes a few hours to get back to you. **Almost all the information you will need to successfully complete the contribution period is enclosed in this document, please read it all before asking questions.**

{% hint style="warning" %}
**Transparency note:** Please note that this internship requires strong knowledge of Python programming language. If you are not yet an expert programmer, we recommend you look for other projects to contribute where your expertise might be better suited.
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

The contribution period runs from <mark style="color:blue;">October 1st</mark> to <mark style="color:blue;">October 29th</mark>. During this time, interested applicants are welcome to contribute to Ersilia's project following the guidelines in this document.

The contribution period is organised in 4 weeks. Each week has a set of specific goals defined, with the objective that mentors can evaluate the intern's experience, interest in the community and team-building work. Once the week's objectives have been met, please focus on:

* Improving your contribution (there is always more publications to read, better bug reports to be written etc)
* Helping out other contributors (we really value group work)

We will be using GitHub issues to track the work of each contributor.

### 📆 WEEK 1: Get to know the community

The first week is focused on getting to know the Ersilia community, our mission and how we work to achieve it. We will be having tons of interactions during the contribution period, so the best is to get to know your peers, your mentors and a bit more about Ersilia.

#### Task 1: Join in the communication channels

**Slack:** we use Slack as our main communication platform, both for the contribution period and afterwards to work with the selected interns. If you have never used this tool, don't worry, is quite intuitive!

1. [Sign up](https://join.slack.com/t/ersilia-outreachy-s24/shared_invite/zt-2e1h1cjs5-lGyrZPprqlpL~qEWgeoCWw) using your preferred email and name.
2. Introduce yourself in the #introductions channel. The general channel is for announcements,  interactions with other fellow contributors, writing tips and suggestions.
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

**Community call:** we will hold a community call on **October 7th 14:00 am CET** to go over the contribution period tasks and answer any questions you might have! The link will be shared via Slack. Attendance is not compulsory, we have tried to find a time that is acceptable for most time zones, we apologize in advance if it means an early start or late end of your day.

**Code of Conduct:** Ersilia is adhered to the [Contributor Covenant Code of Conduct](../../about-us/code-of-conduct.md). Any breaches of the code of conduct, specially harassment or lack of respect for fellow contributors, will mean disqualification as an applicant.

#### Task 2: Install the Ersilia Model Hub

We will be using the Ersilia Model Hub throughout the internship. Please follow the installation [instructions](broken-reference). If you have a UNIX machine (Linux or MacOS) you can install Ersilia directly. If you are using a windows machine you will need a Virtual Machine or a Windows Subsystem Linux (WSL).

{% hint style="info" %}
For Windows users, we recommend using a WSL with Visual Studio Code to access it.
{% endhint %}

{% hint style="danger" %}
A common mistake is to forget the installation of Git-LFS, which is required for many models. Please do so! We also prioritize working with Dockerised models, so make sure to install Docker and Docker Desktop.
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

#### Task 3: Motivation statement

Write, in a thread in your issue, your motivation for joining Outreachy and, in particular, why are you interested in working at Ersilia. A good motivation letter will explain your current skills that are relevant to Ersilia, your reasons to work in Ersilia's project, how this would advance your career and what are your plans during and after the internship.

This needs to be sent to Dhanshree via a private message on Slack

#### Task 4: Obtain approval of the introductory tasks to continue contributing.

If what you have seen and learnt in this first week is appealing and you want to continue working with us, please send a private message to Dhanshree (not any of the other mentors) on Slack. In it:

* Detail which operating system are you using
* Describe your tests of Ersilia and Docker - Demonstrate that you are able to run models and get predictions.&#x20;
* Send the motivation statement to work at Ersilia

{% hint style="warning" %}
We will not be assigning issues or reviewing code contributions of those applicants who have not received the OK to continue working on their application after completing the first week's assignments
{% endhint %}

### 📆 WEEKS 2 and 3: Contribute to a good-first-issue

Once you have successfully completed all the entry-level tasks and received the OK from your mentors to continue contributing, go ahead to  the main repository for the Ersilia Model Hub and have a look at issues labeled as good-first-issue. Please pick an issue where you feel your skills will be best applied. If an issue already has an assigne, you may request to work on it, but only on occasion and depending on the issue we will assign more than one person to the same issue.

<mark style="color:red;">Essential guidelines to work on open issues:</mark>

* Respect the mentors assigning of issues. We will not merge or review any code made by an applicant who was not first assigned to the issue.&#x20;
* Work on your own fork and open a PR when the code is complete and ready for submission.
* Take time to implement the mentors feedback. We might not end up merging the contribution if the quality is not enough.
* Document all the code you write. No documentation will mean no merging of the contribution.

We also recommend looking at other repositories associated to the Ersilia Model Hub for open good issues, such as [ersilia-pack](https://github.com/ersilia-os/ersilia-pack) or [ersilia-self-service](https://github.com/ersilia-os/ersilia-self-service).

### 📆 WEEK 4: Submit your Final Application

Focus the last week <mark style="color:red;">ONLY</mark> in writing your final application to Outreachy. Mentors will not revise any contribution for the last week, only final applications, to ensure we can provide feedback on them.

{% hint style="warning" %}
Final applications must be submitted through the Outreachy website on time. We will not be able to provide help or support for last minute internet connection problems, late submissions and other issues. Please make sure you fill it in with time.
{% endhint %}
