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

* **Phase 1:** Chemistry centered models.
* **Phase 2:** Protein-centered models.
* **Phase 3:** Protein-chemical interaction models.
* **Phase 4:** Cell- and pathogen-based models.
* **Phase 5:** Clinical and epidemiology models.

Within each face, we find models of three types:

* **Type A:** Models developed by others. ****&#x20;
* **Type B:** Models developed by Ersilia based on publicly available data.
* **Type C:** Models developed by Ersilia based on data from collaborators.

{% hint style="info" %}
Ersilia does not follow this roadmap strictly. Since we are already involved in [several projects and collaborations](https://ersilia.io/projects), the Ersilia Model Hub grows in an organic way, necessarily. For example, while we are mainly focused on Phase 1, Type A models at the moment, we spend a lot of time on Phase 1, Type C models, for example as a result of our collaboration with H3D at the University of Cape Town.
{% endhint %}

### Types of models

In other words, we can categorize models

* By input:
  * XX
  * XX
  * XX
* By source
  * XX
  * XX
  * XX

We estimate that, only in Phase 1.A, we will be able to include at least 200 SOTA models. At the end of Phase 1.A, we will write a scientific publication and submit it for peer review.

### Where to find models?

* Scientific literature ([Google Scholar](https://scholar.google.com))
* Code repositories (mainly [GitHub](https://github.com))
* Other model hubs such as [Hugging Face](https://huggingface.co), [TensorFlow Hub](https://tensorflow.org/hub), [PyTorch Hub](https://pytorch.org/hub), etc.
* [Papers With Code](https://paperswithcode.com/)
* Benchmarks, competitions and leaderboards: [Kaggle](https://www.kaggle.com/), [DREAM Challenges](https://dreamchallenges.org/), [Therapeutic Data Commons](https://tdcommons.ai/), etc.
* Twitter: [@KevinKaichuang](https://twitter.com/KevinKaichuang), [@andrewwhite01](https://twitter.com/andrewwhite01)
* Ersilia Slack workspace: the `#literature` channel
* Ersilia models backlog spreadsheet

### Topic of the week

Every **Monday by 2pm CET**, one person of the team will suggest one topic and the rest of us will prioritize model search about this topic during the week. Topics can be related to a task of biomedical relevance, or to a particular family of algorithms. Valid topics could be:

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
Please note that topics of the week should be taken as a **soft guideline**. Discovery and selection of model **should not be blocked** by the existence of the topic of the week. If you find a model that is interesting but it is unrelated to the topic of the week, feel free to select it and work on it. It is still a valid model choice! The topic of the week is just here to help you in case you are lost in your research.
{% endhint %}

#### Announcing the topic of the week

In the Slack `#internships` channel, one person (for example **@Miquel**), will write a message like this:

****:calendar: **@channel** This is the topic of the week!\
:robot: Topic: Antimalarial activity prediction\
:thinking: Why? We are currently working with [Medicines for Malaria Venture](https://mmv.org) and they have asked for predictions.\
:track\_next: Next: **@Gemma**

In this case, **@Miquel** should pin :pushpin: the message so that everybody can find it easily during the week.

It would be great if the rest of the **@channel** can make comments, ask questions, or give feedback about the topic choice. Or even just confirm that you've read the message (:raised\_hands:,:thumbsup:,:heart:,...)

{% hint style="warning" %}
Note that **@Miquel** has nominated **@Gemma**. So **@Gemma** will be responsible for selecting a model in the next week. If you are eager to suggest a topic, simply contact the current responsible present so that they can nominate you :wink:.
{% endhint %}

{% hint style="info" %}
The :thinking: **Why?** bullet point should be short, and can be anything, really. "_I haven't found any model of this kind in the hub"_ is a perfectly valid statement, as is _"I want to learn about this type of models",_ or _"I've read somewhere that there is no drug against this pathogen"_...
{% endhint %}

### How to keep track of models?

The lifecycle of a selected model will follow these three steps:

#### Ersilia Slack `#literature`channel

First, write a quick note in the **Ersilia Slack** `#literature` channel. Simply copy the link to the publication as soon as you discover it, or even a link to a Tweet. Before the link, add the :robot: emoji so that we know it is about a model. For example: :robot: Compound price prediction with deep learning! [\[link\]](https://chemrxiv.org/engage/chemrxiv/article-details/621cf4bace899be245a72621)

#### Ersilia Model Hub Spreadsheet

Then, read about the model in more detail and try to figure out if code is available for it. Then, add the model in a new row in the [**Ersilia Model Hub Spreadsheet**](https://docs.google.com/spreadsheets/d/1TQdei8kkF6zMGyDn0km0qmjZb6p-PM9gsBnSWg3637s/edit?usp=sharing). Please **request edit rights** to **@Gemma** if you don't have them.

In the Spreadsheet, you will find many field. No need to fill all of them! But make sure you provide enough information for others to understand what the model is about and have access to relevant links.

Please write _To do_ in the _Status_ column. Try to provide a _Title_ and a _Description_, and perhaps suggest a _Slug._

#### Ersilia Model Hub AirTable

As soon as we agree to work on a model, we will copy its information from the Spreadsheet to an [AirTable Base](https://airtable.com/shrUcrUnd7jB9ChZV). This Base is accessed programmatically by the [Ersilia CLI](https://github.com/ersilia-os) and our (provisional) [hub interface](https://ersilia.io/model-hub).

{% hint style="success" %}
**You don't have to worry** about the AirTable base for now. As soon as a model is advanced enough, **@Miquel** will add it to the AirTable base.
{% endhint %}

## Protocol



