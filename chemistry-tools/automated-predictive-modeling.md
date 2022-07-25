---
description: Ersilia has developed two AutoML tools that take small molecules as input
---

# Automated predictive modeling

{% hint style="danger" %}
This page is :construction\_worker:**work in progress** :construction\_worker:!
{% endhint %}

### Automated model training

Key to the success of our enterprise is the **quality** of the AI/ML models developed in-house and in collaboration (points 3 and 4 above). To this aim, we plan to use **AutoML** technologies (e.g. hyperparameter optimization), combined with **in-house feature engineering tools**. In brief, our methodology will build upon an extended version of the [Chemical Checker](https://bioactivitysignatures.org) (CC). The CC encapsulates an unprecedented amount of small molecule data in the form of numerical vectors that can be plugged into any standard AI/ML algorithm. In practice, it offers an ideal scenario for **transfer learning**, since it maps small molecules in their relevant bioactivity space, and only a relatively simple fine-tuning or supervised learning step is necessary to provide state-of-the-art predictive tools. Although the CC was initially designed to deal with human cell line data (especially in the context of cancer and Alzheimer’s research), it is easily extensible to **antimicrobial** and **antiviral** data points. We are currently working in this direction with our collaborators.

## Fast predictive modeling with LazyQSAR

### Installation

```bash
# create a conda environment and activate it
conda create -n lazyqsar python=3.7
conda activate lazyqsar
# download lazy-qsar and install it with pip
git clone https://github.com/ersilia-os/lazy-qsar.git
cd lazy-qsar
python -m pip install .
```

### Quick start



## Accurate models with ZairaChem

### Installation

