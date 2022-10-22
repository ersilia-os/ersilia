---
description: Introduction to the Ersilia Model Hub!
---

# Introduction

We are creating the Ersilia Model Hub, a free, **open-source** repository of Artificial Intelligence and Machine Learning (**AI/ML**) models for **drug discovery**. Our platform is aimed at helping researchers identify drug candidates for orphan and neglected diseases, design molecules _de novo_, understand mechanisms of action or anticipate adverse side effects. The ultimate goal of Ersilia is to lower the barrier to drug discovery, encouraging academic groups and companies to pursue the development of new medicines following the principles of **Open Science**. Ersilia disseminates AI/ML models existing in the literature, as well as an in-house collection of models focused on **diseases** that are currently **neglected** by the pharmaceutical industry due to estimated low return on investment. Endemic diseases of low- and middle-income countries (LMIC) belong to this category.

## Technology

You can [browse the Ersilia Model Hub](https://ersilia.io/model-hub) in search of your AI/ML models of interest. All models available in the catalog are accessible through the [Ersilia Python package](https://github.com/ersilia-os/ersilia). At the moment, Ersilia is largely centered on **small molecules**, although we have plans to expand the repository to other areas of biomedical sciences, including genomics, proteomics and epidemiology.

### **Types of models**

The Ersilia Model Hub contains the following AI/ML model categories, sorted by increasing complexity:

1. **Pre-trained literature models**. These correspond to models developed **by others** and published along with code and pre-trained parameters. In this case, our team simply downloads the model and deploys it within the Ersilia Model Hub environment, typically through a Docker container. Unfortunately, according to a preliminary assessment, only a minority of the models fall in this category. A good example is the broad-spectrum antibiotic activity predictor presented in [Stokes et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32084340/).
2. **Re-trained literature models**. In this case, code and data to reproduce the publication results are available, but **pre-trained parameters are not**. This means that re-training is necessary. A substantial number of models are related to this category, including most of the papers testing their models on the [MoleculeNet benchmark](https://pubs.rsc.org/en/content/articlelanding/2018/sc/c7sc02664a).
3. **In-house models trained on relevant datasets**. We are performing an exhaustive literature search effort to **identify drug screening campaigns** relevant to infectious and neglected diseases, like the one reported in [Antonova-Koch et al. (2018)](https://www.science.org/doi/10.1126/science.aat9446), related to a massive screening campaign against the blood-stage malaria parasite. With this data, we develop AI/ML models _de novo_, thanks to our AutoML protocol.
4. **Models built in collaboration**. We partner with experimental scientists who are aligned with our mission and, together, design and develop models **tailored** to answer their specific research questions. One example would be the antimalarial activity predictor that we are developing in the context of the [Open Source Malaria](http://opensourcemalaria.org) project. This predictor is trained on an **experimental dataset provided by the collaborator**, and we are developing a reinforcement learning tool to generate new antimalarial drug candidates and visualize them in the chemical space.

All models developed will be available in our [GitHub repositories](https://github.com/ersilia-os/) under a GPLv3 open source license. Models developed by third parties will include proper **attribution to the original authors**, a link to the source code and a license notice according to the original release.

The figure below provides an overview of the Ersilia Model Hub:

![A schematic view of the Ersilia Model Hub.](https://files.gitbook.com/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-Mj44wxA7bU1hQH19m8I%2Fuploads%2Fgit-blob-6898bac6a4d07a53d774164a0235b77496dfe8a1%2FErsilia\_Hub-02.png?alt=media)

## Roadmap

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

### 1.A: Chemistry models trained by others

Drug discovery, and especially antimicrobial drug discovery, relies heavily on small molecule compounds. **Type A** models can play an important role in the drug discovery pipeline by helping increase the hit rate of pre-clinical experiments, anticipating toxicity outcomes in clinical trials, or even suggesting new drug candidates _de novo_.

In general, the **input** of a Type A model will be a molecules or a list of molecules expressed in the [SMILES](https://en.wikipedia.org/wiki/Simplified\_molecular-input\_line-entry\_system) format. The SMILES strings captures the atomic connectivity of a molecule (i.e. its 2D structure), which is usually sufficient to make good predictions about the compound.

On the contrary, the **output** format can be very variable. In the case of point predictions (for example, [antibiotic activity against _E. coli_](https://github.com/ersilia-os/eos4e40)), the output is one float number (e.g. the IC50). Often, we find multi-output predictions (for example, against the [Tox21 toxicity panel](https://github.com/ersilia-os/eos69p9), containing 12 toxicity outcomes). In this case, the output contains multiple float numbers, each corresponding to one of the outcomes. Moreover, the output can be different than an activity value (a float). For example, in a [generative model](https://github.com/ersilia-os/chem-sampler) the output is one or several molecules, expressed as SMILES, that are generated as derivatives of a seed (input) molecule.

{% hint style="info" %}
If you are only getting started, we highly recommend that you focus on single output float predictions. These typically correspond to compound activities against a certain parasite or protein target of interest.
{% endhint %}
