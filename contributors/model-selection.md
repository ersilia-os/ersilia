---
description: >-
  This page provides guidelines on how to select models to be added in the
  Ersilia Model Hub
---

# Model selection

## Selecting models for the Ersilia Model Hub

### Roadmap

The goal of Ersilia is to create the reference resource of AI/ML models for biomedical research, with a focus on drug discovery for infectious and neglected tropical diseases that affect the Global South.

We have planned the growth of the [Ersilia Model Hub](https://ersilia.io/model-hub) as follows:

* **Phase 1:** Chemistry-centered models.
* **Phase 2:** Protein-centered models.
* **Phase 3:** Protein-chemical interaction models.
* **Phase 4:** Cell- and pathogen-based models.
* **Phase 5:** Clinical and epidemiology models.

Within each face, we find models of three types:

* **Type A:** Models developed by others. ****&#x20;
* **Type B:** Models developed by Ersilia based on publicly available data.
* **Type C:** Models developed by Ersilia based on data from collaborators.

{% hint style="info" %}
Ersilia does not follow this roadmap strictly. Since we are already involved in [several projects and collaborations](https://ersilia.io/projects), the Ersilia Model Hub grows in an organic way, necessarily. For example, while we are currently mainly focused on Phase 1 Type A (**1.A**) models, we spend a lot of time on Phase 1 Type C (**1.C**) models, as a result of our collaboration with key partners like [H3D](http://www.h3d.uct.ac.za/) at the University of Cape Town.
{% endhint %}

We estimate that we will be able to include at least 200 AI/ML models in **1.A**. At the end of this process, we will write a **scientific publication**.

### Phase 1 Type A models: Chemistry models trained by others

The current document focuses on **1.A**, i.e. models developed by third parties having chemical information as input.

Drug discovery, and especially antimicrobial drug discovery, relies heavily on small molecule compounds. Type A models can play an important role in the drug discovery pipeline by helping increase the hit rate of pre-clinical experiments, anticipating toxicity outcomes in clinical trials, or even suggesting new drug candidates _de novo_.

In general, the **input** of a Type A model will be a molecules or a list of molecules expressed in the [SMILES](https://en.wikipedia.org/wiki/Simplified\_molecular-input\_line-entry\_system) format. The SMILES strings captures the atomic connectivity of a molecule (i.e. its 2D structure), which is usually sufficient to make good predictions about the compound.

On the contrary, the **output** format can be very variable. In the case of point predictions (for example, [antibiotic activity against _E. coli_](https://github.com/ersilia-os/eos4e40)), the output is one float number (e.g. the IC50). Often, we find multi-output predictions (for example, against the [Tox21 toxicity panel](https://github.com/ersilia-os/eos69p9), containing 12 toxicity outcomes). In this case, the output contains multiple float numbers, each corresponding to one of the outcomes. Moreover, the output can be different than an activity value (a float). For example, in a [generative model](https://github.com/ersilia-os/chem-sampler) the output is one or several molecules, expressed as SMILES, that are generated as derivatives of a seed (input) molecule.

{% hint style="info" %}
If you are only getting started, we highly recommend that you focus on single output float predictions. These typically correspond to compound activities against a certain parasite or protein target of interest.
{% endhint %}

We will start by populating the Ersilia Model Hub with **1.A** models. We will do this as a team, so that we can all learn from each other. Below, we suggest a procedure to select models and keep track of them.

### Know where to find models and get inspiration

Unfortunately, there is no centralized resource to find AI/ML models for drug discovery online (reason why Ersilia exists). The following resources can help you, though:

* Scientific literature: use [Google Scholar](https://scholar.google.com) and type relevant keywords (e.g. _antimicrobial drug discovery_).
* Reputable scientific journals: identify scientific journals that tend to publish AI/ML research related to drug discovery (e.g. [Chemical Science](https://www.rsc.org/journals-books-databases/about-journals/chemical-science/), [Nature Machine Intelligence](https://www.nature.com/natmachintell/), etc.). Subscribe to their newsletter if you are interested.
* Code repositories (mainly [GitHub](https://github.com)): on a regular Google search, type _github_ next to your keywords for the search.
* Browse other model hubs such as [Hugging Face](https://huggingface.co), [TensorFlow Hub](https://tensorflow.org/hub), [PyTorch Hub](https://pytorch.org/hub), etc.
* Expore [Papers With Code](https://paperswithcode.com/). It is a fantastic resource.
* Follow benchmarks, competitions and leaderboards: [Kaggle](https://www.kaggle.com/), [DREAM Challenges](https://dreamchallenges.org/), [Therapeutic Data Commons](https://tdcommons.ai/), etc.
* Use Twitter: follow accounts like [@KevinKaichuang](https://twitter.com/KevinKaichuang), [@andrewwhite01](https://twitter.com/andrewwhite01), [@EricTopol](https://twitter.com/erictopol) and [@janjensen](https://twitter.com/janjensen). From Ersilia ([@ersiliaio](https://twitter.com/ersiliaio), [@mduranfrigola](https://twitter.com/mduranfrigola), [@TuronGemma](https://twitter.com/TuronGemma)) we tend to retweet relevant literature.
* Use the `#literature` channel in the Ersilia Slack workspace to get inspiration.
* Check the models backlog in the [Ersilia Model Hub Spreadsheet](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=sharing) to avoid duplications and get inspiration too.

### Follow the Topic of the Week

Every **Monday by 9 pm CET**, one person of the team will suggest one topic and the rest of us will prioritize model search about this topic during the week. Topics can be related to a task of biomedical relevance, or to a particular family of algorithms. Valid topics could be:

By biomedical relevance:

* Antimalarial activity prediction
* Broad spectrum antibiotic activity prediction
* Drug toxicity prediction
* Synthetic accessibility of compounds

By algorithmic relevance:

* Graph neural networks
* Reinforcement learning methods
* Chemistry language transformer models

{% hint style="danger" %}
Please note that Topics of the Week should be taken as a **soft guideline**. Discovery and selection of model **should not be blocked** by the existence of the Topic of the Week. If you find a model that is interesting but it is unrelated to the Topic of the Week, feel free to select it and work on it. It is still a valid choice!
{% endhint %}

#### Announcing the Topic of the Week

In the Slack `#internships` channel, one person (for example **@Miquel**) will write a message like this:

****:calendar: **@channel** This is the topic of the week!\
:robot: Topic: Antimalarial activity prediction\
:thinking: Why? We are currently working with [Medicines for Malaria Venture](https://mmv.org) and they have asked for predictions on some antimalarial candidates.\
:track\_next: Next: **@Gemma**

In this case, **@Miquel** should pin :pushpin: the message so that everybody can find it easily during the week.

It would be great if the rest of the **@channel** can make comments, ask questions, or give feedback about the topic choice. Or even just confirm that you've read the message (:raised\_hands:,:thumbsup:,:heart:,...)

{% hint style="warning" %}
Note that **@Miquel** has nominated **@Gemma**. So **@Gemma** will be responsible for selecting a model next week. If you are eager to suggest a topic, simply contact the current responsible person so that they can nominate you :wink:.
{% endhint %}

{% hint style="info" %}
The :thinking: **Why?** bullet point should be short, and can be anything, really. "_I haven't found any model of this kind in the hub"_ is a perfectly valid statement, as is _"I want to learn about this family of models",_ or _"I've read somewhere that there is no drug against this pathogen"_, etc.
{% endhint %}

### Notify and keep track of models

It is important that you communicate your research to the rest of the team. We suggest the following three steps:

#### Ersilia Slack `#literature`channel

First, write a quick note in the **Ersilia Slack** `#literature` channel. Simply copy the link to the publication as soon as you discover the model, or even a link to a tweet. Before the link, add the :robot: emoji so that we know it is about a model. For example: :robot: _Compound price prediction with deep learning!_ [_\[link\]_](https://chemrxiv.org/engage/chemrxiv/article-details/621cf4bace899be245a72621)__

{% hint style="info" %}
The Slack `#literature` channel contains models :robot: and much more. We encourage you to check this channel regularly and be active in it! This is the best way to share our collective knowledge and avoid duplication of efforts. Naturally, this channel will not be very structured, and we feel it shouldn't be. Rather, it should be a space to quickly share findings and find inspiration for your research.
{% endhint %}

#### Ersilia Model Hub Spreadsheet

Once you've posted the model in the `#literature` channel, read about it in more detail and try to figure out if code is available. Then, add the model as a new entry in the [**Ersilia Model Hub Spreadsheet**](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=sharing). Please **request edit rights** to **@Gemma** if you don't have them.

In the Spreadsheet, you will find two sheets:

* **Hub:** Contains complete documentation about the models, including an Ersilia Open Source identifier (e.g. `eos4e40`), a slug (e.g. `chemprop-antibiotic`) and a status (e.g. _Done_).
* **Raw List:** Contains a backlog of models that could be of interest to Ersilia. This sheet contains minimal information about the models. There is an _Approved_ and a _Selected_ tickbox. _Approved_ means that Ersilia is willing to incorporate this model. _Selected_ means that someone is already working on the model or the model has been successfully incorporated in the Ersilia Model Hub.

You should **start by the Raw List sheet**. You are always free to add models there (relevant to the [Topic of the Week](model-selection.md#follow-the-topic-of-the-week) or not). **@Miquel** is responsible for curating this list and approving the models, he will use the Slack `#internships` channel to ask questions or discuss the relevance of the model before approving them.

Start by adding your model in the Raw List sheet and wait for approval. As soon as the model has been approved, you are ready to move forward. **Please tick the **_**Selected**_** column** so that nobody else picks your model of interest.

Now you are ready to start **filling the Hub sheet**. Fill in as much information as possible, and if something is missing make sure you provide relevant links and enough information for others to understand what the model is about.

Please write _To do_, _In progress_ or _Done_ in the _Status_ column:

* _To do_ means that you have filled the information but you have not started the model incorporation per se. In other words, you have not started to work on the coding part yet.
* _In progress_ means that you are already working on the code.
* _Done_ means that the model has been successfully incorporated in the Ersilia Model Hub.

Try to provide a _Title_ and a _Description_, and perhaps suggest a _Slug._ The Ersilia Open Source identifiers are already predefined, so no need to worry about it. Don't forget to write your name in the _Contributor_ column.

{% hint style="info" %}
**@Gemma** and **@Ife** are the maintainers of the Spreadsheet. Please reach out to them if you have questions or suggestions.
{% endhint %}

{% hint style="success" %}
The goal of the _Contributor_ colum is, simply, to know who is the person of reference for each model. Some models are more difficult to add than others, so **you should not be stressed** about the number of models that you contribute. Ersilia is a safe and collaborative space, we do not monitor this kind of metrics :hugging:.
{% endhint %}

#### Ersilia Model Hub AirTable

For now, **@Miquel** will be responsible for **curating** the models listed in the Spreadsheet. He will reach out to you if he has questions or needs more information about the model. He will copy the information from the Spreadsheet to an [AirTable Base](https://airtable.com/shrUcrUnd7jB9ChZV). This Base is accessed programmatically by the [Ersilia CLI](https://github.com/ersilia-os) and our (provisional) [hub interface](https://ersilia.io/model-hub).

{% hint style="success" %}
**You don't have to worry** about the AirTable base for now. This database is fully managed by **@Miquel**.
{% endhint %}

### Start coding! :rocket:

As soon as you feel ready to start coding, you should change the _Status_ in the Spreadsheet Hub sheet from _To do_ to _In progress_. Please go to the next page to learn more about how to [Incorporate models to the hub](incorporate-models.md).

{% content-ref url="incorporate-models.md" %}
[incorporate-models.md](incorporate-models.md)
{% endcontent-ref %}

## TL;DR

In brief, this is a suggested routine that you can follow:

1. Check the [Topic of the Week](model-selection.md#follow-the-topic-of-the-week).
2. [Search models](model-selection.md#know-where-to-find-models-and-get-inspiration) related to the topic.
3. [Post your findings](model-selection.md#ersilia-slack-literaturechannel) in the `#literature` channel.
4. Choose one or few models and [add them in the backlog](model-selection.md#ersilia-model-hub-spreadsheet) (**Raw List** sheet of the [Ersilia Model Hub Spreadsheet](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=sharing)).
5. Wait for approval.
6. Select an approved model.
7. Move to the **Hub** sheet of the Spreadsheet and [fill in as much information as you can](model-selection.md#ersilia-model-hub-spreadsheet). For now, set the _Status_ to _To do_.
8. When you are ready to [start coding](incorporate-models.md), change the _Status_ to _In progress_.
9. When the model is successfully added, change the _Status_ to _Done_! __ :tada:__
