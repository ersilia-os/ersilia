---
description: This tutorial explains how to incorporate models in the Ersilia Model Hub
---

# Contribute models

This page serves as a guideline to develop a selected model to add to the Ersilia Model Hub. As [described earlier](../ersilia-model-hub/introduction.md), models can be of three types, namely:

1. Models developed by third parties.
2. Models developed by Ersilia based on publicly available data.
3. Models developed by Ersilia based on data from collaborators.

In this page, we focus on the first type of models, i.e. models developed by **third parties**. These third party models come with source code and parameters (or data to train them). We show how to to use the Ersilia Model Template to adapt the assets provided by third parties and incorporate them into the Ersilia Model Hub.

## Model List

First, check our Model Hub Backend, which contains:

* **Model Hub:** [models](https://airtable.com/shrNc3sTtTA3QeEZu/tblZGe2a2XeBxrEHP) that are already in the Hub or being incorporated at the moment. Each model should contain accurate information about its usage, source and any other relevant information.
* **Model Suggestions:** Contains a [backlog](https://airtable.com/shrZdogkgUk8ICyO5) of models that could be of interest to Ersilia. This sheet contains minimal information about the models. There is an _Approved_ and a _Selected_ tickbox. _Approved_ means that Ersilia is willing to incorporate this model. _Selected_ means that someone is already working on the model or the model has been successfully incorporated in the Ersilia Model Hub.

To start contributing:

1. Become a member of the **Ersilia Slack** Workspace. All our communications are centralised there. If you need a sign up link, please email us at hello@ersilia.io
2. Join the **#contributors** channel in Slack and ask for edit permits for the Spreadsheet.
3. Fill in models of interest in the **Model Suggestions** list (if you need some sources of inspiration, check below)
4. Wait for model to be **approved** (or select an already approved one) and notify everyone of your selection on the #contributors channel to avoid duplication.
5. Fill in the information on the Model Hub sheet. Use the next available EOS ID and **do not** modify it.
6. Start coding! Head to the [Model Incorporation Guidelines](../ersilia-model-hub/contribute-models/model-incorporation-guidelines.md) for more information.

### Know where to find models and get inspiration

Unfortunately, there is no centralized resource to find AI/ML models for drug discovery online (reason why Ersilia exists). The following resources can help you, though:

* Scientific literature: use [Google Scholar](https://scholar.google.com) and type relevant keywords (e.g. _antimicrobial drug discovery_).
* Reputable scientific journals: identify scientific journals that tend to publish AI/ML research related to drug discovery (e.g. [Chemical Science](https://www.rsc.org/journals-books-databases/about-journals/chemical-science/), [Nature Machine Intelligence](https://www.nature.com/natmachintell/), etc.). Subscribe to their newsletter if you are interested.
* Code repositories (mainly [GitHub](https://github.com)): on a regular Google search, type _github_ next to your keywords for the search.
* Browse other model hubs such as [Hugging Face](https://huggingface.co), [TensorFlow Hub](https://tensorflow.org/hub), [PyTorch Hub](https://pytorch.org/hub), etc.
* Expore [Papers With Code](https://paperswithcode.com/). It is a fantastic resource.
* Follow benchmarks, competitions and leaderboards: [Kaggle](https://www.kaggle.com/), [DREAM Challenges](https://dreamchallenges.org/), [Therapeutic Data Commons](https://tdcommons.ai/), etc.
* Use Twitter: follow accounts like [@KevinKaichuang](https://twitter.com/KevinKaichuang), [@andrewwhite01](https://twitter.com/andrewwhite01), [@EricTopol](https://twitter.com/erictopol), [@iamtrask](https://twitter.com/iamtrask) and [@janjensen](https://twitter.com/janjensen). From Ersilia ([@ersiliaio](https://twitter.com/ersiliaio), [@mduranfrigola](https://twitter.com/mduranfrigola), [@TuronGemma](https://twitter.com/TuronGemma)) we tend to retweet relevant literature.
* Use the `#literature` channel in the Ersilia Slack workspace to get inspiration.
* Check the models backlog in the [Ersilia Model Hub](https://airtable.com/shrNc3sTtTA3QeEZu) to avoid duplications and get inspiration too.


