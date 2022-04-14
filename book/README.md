# 🤗 Welcome to Ersilia!

## The Ersilia Model Hub

We are creating the Ersilia Model Hub, a free, **open-source** repository of Artificial Intelligence and Machine Learning (**AI/ML**) models for **drug discovery**. Our platform is aimed at helping researchers identify drug candidates for orphan and neglected diseases, design molecules _de novo_, understand mechanisms of action or anticipate adverse side effects. The ultimate goal of Ersilia is to lower the barrier to drug discovery, encouraging academic groups and companies to pursue the development of new medicines following the principles of **Open Science**. Ersilia disseminates AI/ML models existing in the literature, as well as an in-house collection of models focused on **diseases** that are currently **neglected** by the pharmaceutical industry due to estimated low return on investment. Endemic diseases of low- and middle-income countries (LMIC) belong to this category.

## Technology

You can [browse the Ersilia Model Hub](https://airtable.com/shrUcrUnd7jB9ChZV) in search of your AI/ML models of interest. All models available in the catalog are accessible through the [Ersilia Python package](https://github.com/ersilia-os/ersilia). At the moment, Ersilia is largely centered on **small molecules**, although we have plans to expand the repository to areas of biomedicine, including genomics and epidemiology.

### **Types of models**

The Ersilia Model Hub contains the following AI/ML model categories, sorted by increasing complexity:

1. **Pre-trained literature models**. These correspond to models developed **by others** and published along with code and pre-trained parameters. In this case, our team simply downloads the model and deploys it within the Ersilia Model Hub environment, typically through a Docker container. Unfortunately, according to a preliminary assessment, only a minority of the models fall in this category. A good example is the broad-spectrum antibiotic activity predictor presented in [Stokes et al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32084340/).
2. **Re-trained literature models**. **** In this case, code and data to reproduce the publication results are available, but **pre-trained parameters are not**. This means that re-training is necessary. A substantial number of models are related to this category, including most of the papers testing their models on the [MoleculeNet benchmark](https://pubs.rsc.org/en/content/articlelanding/2018/sc/c7sc02664a).
3. **In-house models trained on relevant datasets**. We are performing an exhaustive literature search effort to **identify drug screening campaigns** relevant to infectious and neglected diseases, like the one reported in [Antonova-Koch et al. (2018)](https://www.science.org/doi/10.1126/science.aat9446), related to a massive screening campaign against the blood-stage malaria parasite. With this data, we develop AI/ML models _de novo_, thanks to our AutoML protocol.
4. **Models built in collaboration**. We partner with experimental scientists who are aligned with our mission and, together, design and develop models **tailored** to answer their specific research questions. One example would be the antimalarial activity predictor that we are developing in the context of the [Open Source Malaria](http://opensourcemalaria.org) project. This predictor is trained on an **experimental dataset provided by the collaborator**, and we are developing a reinforcement learning tool to generate new antimalarial drug candidates and visualize them in the chemical space.

All models developed will be available in our [GitHub repositories](https://github.com/ersilia-os/) under a permissive open source license, typically MIT or CC-BY-4.0. Models developed by third parties will include proper **attribution to the original authors**, a link to the source code and a license notice according to the original release.

### Automated model training

Key to the success of our enterprise is the **quality** of the AI/ML models developed in-house and in collaboration (points 3 and 4 above). To this aim, we plan to use **AutoML** technologies (e.g. hyperparameter optimization), combined with **in-house feature engineering tools**. In brief, our methodology will build upon an extended version of the [Chemical Checker](https://bioactivitysignatures.org) (CC). The CC encapsulates an unprecedented amount of small molecule data in the form of numerical vectors that can be plugged into any standard AI/ML algorithm. In practice, it offers an ideal scenario for **transfer learning**, since it maps small molecules in their relevant bioactivity space, and only a relatively simple fine-tuning or supervised learning step is necessary to provide state-of-the-art predictive tools. Although the CC was initially designed to deal with human cell line data (especially in the context of cancer and Alzheimer’s research), it is easily extensible to **antimicrobial** and **antiviral** data points. We are currently working in this direction with our collaborators.

The figure below provides an overview of the Ersilia Model Hub:

![A schematic view of the Ersilia Model Hub.](<.gitbook/assets/Ersilia\_Hub-02 (1).png>)

* To explore our **open source contributions**, visit our [GitHub profile](https://github.com/ersilia-os/).
* To learn more **about our organisation**, visit the [Ersilia Open Source Initiative](https://ersilia.io) website.
* If you want to know more or contribute, please **reach out** to us at [hello@ersilia.io](mailto:hello@ersilia.io)!
