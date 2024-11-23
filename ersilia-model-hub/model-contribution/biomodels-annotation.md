---
description: >-
  This page describes how to annotate Ersilia models in the BioModels Tool
  contributing towards FAIRness.
---

# BioModels Annotation

## Background

Sharing of machine learning models most importantly in the field of drug discovery is important in creating a FAIReR (Findable, Accessible, Interoperable, Reusable, and Reproducible) collection of machine learning models which in turn makes it easier to reproduce and reuse these models. This reduces the need to rebuild models from scratch, increases their usefulness in various applications, and speeds up progress in drug discovery.

In the process of making ML models FAIR and shareable, there are standards and protocols to follow which includes; sharing model training code, dataset information, reproduced figures, model evaluation metrics, trained models, Docker files, model metadata, and FAIR dissemination. Here, we will focus on the Model Metadata, and its annotation.

## Model Metadata

To share ML models effectively, it’s important to provide relevant information about the models. The information being the Metadata. Metadata is organised information that describes, explains, or helps find, use, or manage a resource. In the context of Ersilia Models, metadata is data about the model and it is classified into three categories namely; **Biological Metadata**, **Computational Metadata**, and the **Description Metadata**. The metadata enables the findability and accessibility of the models based on its specific characteristics by other researchers and modellers.

*   #### Biological Metadata

    The biological relevance of a model is an important aspect of a model. This ranges from its bioactivity, biological processes explained by the model, biological system the model was trained on, tissue or cell type involved, assay type, biological entity, Ersilia model theme (ranges from infectious disease to ADME property) and Compartment in which biological process is happening.
*   #### Computational Metadata

    The metadata identifies the model based on its specific characteristics, such as the type of ML algorithm used, the modelling approach, its evaluation metrics, and the functional properties of the model such as input data type, and model output. &#x20;
*   #### Description Metadata

    A model is described by its publication, a code base such as GitHub, data repository such as Zenodo, and lastly its deployment which could be in the form of a web server.&#x20;

## Why is Annotation Important?

An annotation is an association of a metadata with an ontology term. Annotating a model consists of mapping the identified model metadata with terms from controlled vocabularies and entries in data resources.&#x20;

Annotation of a Model Metadata are crucial to:

* Precisely identify model categories
* Improve understanding of the model's structure
* Make it easier to compare different models
* Simplify model integration
* Enable efficient searches
* Add meaningful context to the model
* Enhance understanding of the biology behind the model
* Allow conversion and reuse of the model
* Facilitate integration of the model with biological knowledge

### Controlled Vocabularies and Ontology

Controlled vocabularies contain set terms that describe concepts in a specific domain. These terms have definitions that help us understand and agree on their meanings, and also helps with indexing and easy retrieval. Ontologies use controlled vocabularies to describe concepts and how they are related in a structured & computable format.&#x20;

**The following ontologies are preferred**;

* [NCI Thesaurus OBO Edition NCIT](https://www.ebi.ac.uk/ols4/ontologies/ncit)
* [BioAssay Ontology BAO](https://www.ebi.ac.uk/ols4/ontologies/bao)
* [Bioinformatics Concept EDAM](https://www.ebi.ac.uk/ols4/ontologies/edam)
* [STATO: the statistical methods ontology](https://www.ebi.ac.uk/ols4/ontologies/stato)
* [The BRENDA Tissue Ontology (BTO)](https://www.ebi.ac.uk/ols4/ontologies/bto)
* [Chemical Entities of Biological Interest CHEBI](https://www.ebi.ac.uk/ols4/ontologies/chebi)
* [Chemical Information Ontology (cheminf)](https://www.ebi.ac.uk/ols4/ontologies/cheminf)
* [Chemical Methods Ontology CHMO](https://www.ebi.ac.uk/ols4/ontologies/chmo)
* [Drug-drug Interaction and Drug-drug Interaction Evidence Ontology (DIDEO) ](https://www.ebi.ac.uk/ols4/ontologies/dideo)
* [Gene Ontology GO](https://www.ebi.ac.uk/ols4/ontologies/go)
* [Infectious Disease Ontology IDO](https://www.ebi.ac.uk/ols4/ontologies/ido)
* [Model Card Ontology MCRO](https://www.ebi.ac.uk/ols4/ontologies/mcro)
* [Molecular Interactions MI](https://www.ebi.ac.uk/ols4/ontologies/mi)
* [Mondo Disease Ontology MONDO](https://www.ebi.ac.uk/ols4/ontologies/mondo)
* [NCBI Taxonomy](https://www.ebi.ac.uk/ols4/ontologies/ncbitaxon)
* [OBCS: Ontology of Biological and Clinical Statistics](https://www.ebi.ac.uk/ols4/ontologies/obcs)
* [Ontology for Biomedical Investigations OBI](https://www.ebi.ac.uk/ols4/ontologies/obi)
* [Ontology for MIRNA Target OMIT](https://www.ebi.ac.uk/ols4/ontologies/omit)
* [Ontology for Parasite Lifecycle OPL](https://www.ebi.ac.uk/ols4/ontologies/opl)
* [PATO - the Phenotype And Trait Ontology](https://www.ebi.ac.uk/ols4/ontologies/pato)
* [PRotein Ontology PR](https://www.ebi.ac.uk/ols4/ontologies/pr)
* [Uber-anatomy ontology UBERON](https://www.ebi.ac.uk/ols4/ontologies/uberon)
* [Experimental Factor Ontology EFO](https://www.ebi.ac.uk/ols4/ontologies/efo)

## How do we Annotate?

This is the process of curating an annotation file with all essential information. Annotation is done by linking the right ontology, a cross-reference to a metadata with the addition of values and qualifiers in order to explicitly define the relationship between the model metadata, and the linked resources.&#x20;

1. Metadata is the data of a data. The latter being the model. This model is referred to as the entity in annotation.
   * Entity (eg. **Model**)
2. Extract all information available for each metadata categories.
3. Each annotation is linked to external data resources and values e.g., EDAM, STATO. An external data resource could be a database of ontology or the ontology itself.&#x20;
   * This improves the model quality
   * Essential for the model search criteria
4. **Value** enhances the accessibility and integrates a metadata with other data resources using a compact identifier.
   * **Metadata** (eg. Machine learning)
   * **Ontology** (Bioinformatics Concept EDAM: edam:topic\_3474)
   * **Value** (https://identifiers.org/bptl/edam:topic\_3474)
5.  Append qualifiers to each annotation. [**Qualifiers**](https://drive.google.com/file/d/1JqjcH0T0UTWMuBj-scIMwsyt2z38A0vp/view) explain the relationship between a metadata and the model itself.&#x20;

    * **Qualifier** (eg. **bqbiol** and **bqmodel**)
    * A relationship is either biological - bqbiol, or computational - bqmodel.&#x20;

    The following qualifiers are used to describe relationships in the annotation;

    * **bqbiol:hasTaxon** - describes a relationship between a model and organism
    * **bqbiol:occursIn** - a compartment where a process occurs
    * **bqbiol:hasProperty** - general biological property
    * **bqmodel:hasProperty** - all model properties
    * **bqmodel:isDescribedBy** - the model resources
    * **bqbiol:hasInput** - model input data
    * **bqbiol:hasDataset** - the model training data
    * **bqbiol:hasOutput** - biological output of the model
6. [**Dome Annotation**](https://www.nature.com/articles/s41592-021-01205-4) gives more context to the computational metadata&#x20;
   * **D - Data** ( this could be data source or the type of input data
   * **O - Optimization** (each model has their algorithm)
   * **M - Model** (the model source code and it’s executable form)
   * **E - Evaluation metrics** (all model performance are evaluated)&#x20;

## Tools for Annotation

### [Ontology Lookup Service](https://www.ebi.ac.uk/ols4/ontologies)

**OLS** is a search and visualisation service that hosts 260+ biological and biomedical ontologies in one place.&#x20;

### [Zooma](https://www.ebi.ac.uk/spot/zooma/)

**Zooma** is a ontology mapping tool that can be used to automatically map free text&#x20;

## Steps to Annotating a Model - An Example

### 1. Identify a Model and its associated Information

Associated information includes the model publication, repository and its source code.&#x20;

#### Read the Publication

All model metadata are enclosed in its publication and it’s important to read the publication to understand the biological or chemical processes the model performs, its bioactivity, the algorithm of the model, its training data and the validation performed among other information. Go through the repository to validate the information in the publication.

#### Special Scenario

There are models that are built by fine-tuning other large models with different datasets, and performing several tasks. It’s important to understand the base models in this case, and all its properties.

#### Case-study

For the purpose of example, we're working with an [**antimalarial model**](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00487-2) with the tag [**eos4zfy**](https://github.com/ersilia-os/eos4zfy) from the Ersilia Model. In this case, this is a collaborative project between the EMBL-EBI and other big pharma and institutes. Here, an individual QSAR model to identify novel molecules that may have antimalarial properties built on private dataset was merged together to develop MAIP. A free web platform available for mass prediction of potential malaria inhibiting compounds.&#x20;

### 2. Assign the Metadata Entity

The metadata entity is the source of the metadata. It’s more of the Metadata Data which is its Model. Adding the metadata entity makes the table looks like this;

<table><thead><tr><th>Entity</th><th data-hidden></th><th data-hidden></th></tr></thead><tbody><tr><td>Model</td><td></td><td></td></tr><tr><td>Model</td><td></td><td></td></tr><tr><td>Model</td><td></td><td></td></tr></tbody></table>

### 3. Extract all information available for each metadata category.

To identify all model metadata associated with this model. We’d go through the publication, the web platform to understand its pipeline, its source code and Ersilia implementation process. This [**template**](https://docs.google.com/spreadsheets/d/1bsCkN5Ugmo3tSF4oc-NGAQaZYfTP9LIABezrov-Vj-0/edit#gid=0) can be adapted to individual use and a sample visual can be seen below.

#### Biological Metadata

These metadata are extracted from the abstract section of the publication and includes the disease, causative agents, data classification and the biological output of the model.&#x20;

#### Computational Metadata

These includes the algorithm of the model, it’s evaluation method, and data type

#### Descriptive Metadata

Each model is described by a publication, source code and its implementation.&#x20;

| Entity | Model Metadata Categories | Metadata                                                                                                                                                                               |
| ------ | ------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Model  | Biological Metadata       | <ol><li>Homo sapiens</li><li>Plasmodium falciparum</li><li>Malaria</li><li>Antimalarial properties</li><li>Active</li><li>Inactive</li><li>Antimalarial compounds prediction</li></ol> |
| Model  | Computational Metadata    | <ol start="8"><li>Classification models</li><li>Naïve Bayesian model</li><li>AUC–ROC</li><li>5-fold cross validation</li><li>Smiles descriptors</li><li>malaria dataset</li></ol>      |
| Model  | Descriptive Metadata      | <ol start="14"><li>MAIP web platform - Source code</li><li>Ersilia Incorporation URL</li><li>PubMed URL</li></ol>                                                                      |

For the purpose of example, these are sample Metadata from this model and its classification.

_**P.S: Column 2 (Model Metadata Categories) is just for descriptive purpose. It's not part of the annotation**_

#### Special Scenario

Some models were validated after building experimentally. This validation is done either in-vivo or in-vitro. It occurs as a predicted compound from the model being further validated experimentally to confirm its bioactivity. This validation is important for models that undergo such and should be annotated for the model.&#x20;

### 4. Map the Metadata to the right Ontology

This is the main process of annotation, and it’s associating a metadata to the right ontology. Ontology can be identified through the [**Ontology Lookup Service**](https://www.ebi.ac.uk/ols4). The Ontology Lookup Service (OLS) is a repository for biomedical ontologies that aims to provide a single point of access to the latest ontology versions.&#x20;

To ensure standardization and interoperability, it's crucial to identify relevant ontologies through the Ontology Lookup Service (OLS). These ontologies will help in annotating the model components accurately.  Search for your terms, for example, Machine Learning, in the search bar, and select the right term in the preferred ontology. If not found in the preferred ontology, look through other available options with the right meaning.

{% hint style="info" %}
Sometimes, the exact term isn't found in the OLS, and in this case, the closest term can be used to replace the metadata.&#x20;
{% endhint %}

In choosing the right ontology for a metadata, there are important things to consider.&#x20;

1. The ontology with the best metadata meaning&#x20;
2. &#x20;Inclusive of a preferred ontology for better indexing

After mapping the metadata to the right ontology, we have a table like this;

| Entity | Preferred Ontology                                                                                                                                                                                                                                                                                            | Metadata                                                                                                                                                                               |
| ------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Model  | <ol><li>NCBI Taxonomy</li><li>NCBI Taxonomy</li><li>Experimental Factor Ontology EFO</li><li>NCI Thesaurus OBO Edition NCIT</li><li>NCI Thesaurus OBO Edition NCIT</li><li>NCI Thesaurus OBO Edition NCIT</li><li>Chemical Entities of Biological Interest CHEBI</li></ol>                                    | <ol><li>Homo sapiens</li><li>Plasmodium falciparum</li><li>Malaria</li><li>Antimalarial properties</li><li>Active</li><li>Inactive</li><li>Antimalarial compounds prediction</li></ol> |
| Model  | <ol start="8"><li>STATO: the statistical methods ontology</li><li>STATO: the statistical methods ontology</li><li>STATO: the statistical methods ontology</li><li>Ontology for Biomedical Investigations OBI</li><li>Chemical information ontology (cheminf)</li><li>NCI Thesaurus OBO Edition NCIT</li></ol> | <ol start="8"><li>Classification models</li><li>Naïve Bayesian model</li><li>AUC–ROC</li><li>5-fold cross validation</li><li>Smiles descriptors</li><li>malaria dataset</li></ol>      |
| Model  | **Special cases**                                                                                                                                                                                                                                                                                             | <ol start="14"><li>MAIP web platform - Source code</li><li>Ersilia Incorporation URL</li><li>PubMed URL</li></ol>                                                                      |



{% hint style="info" %}
Note that the descriptive Medatada, like the PubMed URL of the model, or the specific Ersilia GitHub repository where the model is hosted, do not have an ontology (special cases).&#x20;
{% endhint %}

### 5. Link the Ontology to their values

Values enhance accessibility and integrate metadata with other data resources in the form of a URL (Uniform Resource Locator).

Each ontology has its **accession identifier** and a value is formed using the ontology identifier with a **compact identifier**. The compact identifier is a resolution service that provides consistent access in form of [https://identifiers.org/](https://identifiers.org/)

1. Value = [https://identifiers.org/](https://identifiers.org/) + NCIT:C176231
2. Value = [https://identifiers.org/NCIT:C176231](https://identifiers.org/NCIT:C176231)&#x20;

Each ontology is linked to their respective value using the formula above, and the table looks like this;

<table><thead><tr><th width="97">Entity</th><th>Preferred Ontology</th><th width="220">Values</th><th>Metadata</th></tr></thead><tbody><tr><td>Model</td><td><ol><li>NCBI Taxonomy</li><li>NCBI Taxonomy</li><li>Experimental Factor Ontology EFO</li><li>NCI Thesaurus OBO Edition NCIT</li><li>NCI Thesaurus OBO Edition NCIT</li><li>NCI Thesaurus OBO Edition NCIT</li><li>Chemical Entities of Biological Interest CHEBI</li></ol></td><td><ol><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068">https://identifiers.org/taxonomy/9606</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/taxonomy/5833</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/EFO:0001068</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/NCIT:C271</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068">https://identifiers.org/NCIT:C45329</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/NCIT:C154407</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/CHEBI:38068</a></li></ol></td><td><ol><li>Homo sapiens</li><li>Plasmodium falciparum</li><li>Malaria</li><li>Antimalarial properties</li><li>Active</li><li>Inactive</li><li>Antimalarial compounds prediction</li></ol></td></tr><tr><td>Model</td><td><ol start="8"><li>STATO: the statistical methods ontology</li><li>STATO: the statistical methods ontology</li><li>STATO: the statistical methods ontology</li><li>Ontology for Biomedical Investigations OBI</li><li>Chemical information ontology (cheminf)</li><li>NCI Thesaurus OBO Edition NCIT</li></ol></td><td><ol start="8"><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/STATO:0000031<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/STATO:0000530<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/STATO:0000274<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/obi:OBI_0200032<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/CHEMINF:000018<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/NCIT:C47824</a></li></ol></td><td><ol start="8"><li>Classification models</li><li>Naïve Bayesian model</li><li>AUC–ROC</li><li>5-fold cross validation</li><li>Smiles descriptors</li><li>malaria dataset</li></ol></td></tr><tr><td>Model</td><td><ol start="14"><li>Online Web server</li><li>Ersilia Model Hub</li><li>PubMed Identification Number PMID</li></ol></td><td><ol start="14"><li><a href="https://www.ebi.ac.uk/chembl/maip/https://github.com/ersilia-os/eos4zfyhttps://identifiers.org/pubmed:33618772">https://www.ebi.ac.uk/chembl/maip/<br></a></li><li><a href="https://www.ebi.ac.uk/chembl/maip/https://github.com/ersilia-os/eos4zfyhttps://identifiers.org/pubmed:33618772">https://github.com/ersilia-os/eos4zfy<br></a></li><li><a href="https://www.ebi.ac.uk/chembl/maip/https://github.com/ersilia-os/eos4zfyhttps://identifiers.org/pubmed:33618772">https://identifiers.org/pubmed:33618772</a></li></ol></td><td><ol start="14"><li>MAIP web platform - Source code</li><li>Ersilia Incorporation URL</li><li>PubMed URL</li></ol></td></tr></tbody></table>

### 6. Associate the right qualifier to each annotation

Each metadata as previously explained is either a biology component of the model or a computational component or a descriptive component. Here, we’d annotate the metadata based on the category the fall.

For example;

| Metadata Category      | Metadata                  | Qualifier             |
| ---------------------- | ------------------------- | --------------------- |
| Biological Metadata    | Malaria                   | bqbiol:hasProperty    |
| Computational Metadata | naïve Bayesian model      | bqmodel:hasProperty   |
| Descriptive Metadata   | Ersilia Incorporation URL | bqmodel:isDescribedBy |

After adding qualifiers to each metadata, the table looks like this;

<table><thead><tr><th width="97">Entity</th><th>Qualifiers</th><th>Preferred Ontology</th><th width="220">Values</th><th>Metadata</th></tr></thead><tbody><tr><td>Model</td><td><ol><li>bqbiol:hasTaxon</li></ol><ol start="2"><li>bqbiol:hasTaxon</li></ol><ol start="3"><li>bqbiol:hasProperty</li></ol><ol start="4"><li>bqbiol:hasProperty</li></ol><ol start="5"><li>bqbiol:hasProperty</li></ol><ol start="6"><li>bqbiol:hasProperty</li></ol><ol start="7"><li>bqbiol:hasOutput</li></ol><p><br></p></td><td><ol><li>NCBI Taxonomy</li><li>NCBI Taxonomy</li><li>Experimental Factor Ontology EFO</li><li>NCI Thesaurus OBO Edition NCIT</li><li>NCI Thesaurus OBO Edition NCIT</li><li>NCI Thesaurus OBO Edition NCIT</li><li>Chemical Entities of Biological Interest CHEBI</li></ol></td><td><ol><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068">https://identifiers.org/taxonomy/9606</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/taxonomy/5833</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/EFO:0001068</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/NCIT:C271</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068">https://identifiers.org/NCIT:C45329</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/NCIT:C154407</a></li><li><a href="https://identifiers.org/taxonomy/9606https:/identifiers.org/taxonomy/5833https:/identifiers.org/EFO:0001068https:/identifiers.org/NCIT:C271https:/identifiers.org/NCIT:C45329https:/identifiers.org/NCIT:C154407https:/identifiers.org/CHEBI:38068"><br>https://identifiers.org/CHEBI:38068</a></li></ol></td><td><ol><li>Homo sapiens</li><li>Plasmodium falciparum</li><li>Malaria</li><li>Antimalarial properties</li><li>Active</li><li>Inactive</li><li>Antimalarial compounds prediction</li></ol></td></tr><tr><td>Model</td><td><ol start="8"><li>bqmodel:hasProperty</li></ol><ol start="9"><li>bqmodel:hasProperty</li></ol><ol start="10"><li>bqmodel:hasProperty</li></ol><ol start="11"><li>bqmodel:hasProperty</li></ol><ol start="12"><li>bqbiol:hasInput</li></ol><ol start="13"><li>bqbiol:hasDataset</li></ol><p><br></p></td><td><ol start="8"><li>STATO: the statistical methods ontology</li><li>STATO: the statistical methods ontology</li><li>STATO: the statistical methods ontology</li><li>Ontology for Biomedical Investigations OBI</li><li>Chemical information ontology (cheminf)</li><li>NCI Thesaurus OBO Edition NCIT</li></ol></td><td><ol start="8"><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/STATO:0000031<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/STATO:0000530<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/STATO:0000274<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/obi:OBI_0200032<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/CHEMINF:000018<br></a></li><li><a href="https://identifiers.org/STATO:0000031https:/identifiers.org/STATO:0000530https:/identifiers.org/STATO:0000274https:/identifiers.org/obi:OBI_0200032https:/identifiers.org/CHEMINF:000018https:/identifiers.org/NCIT:C47824">https://identifiers.org/NCIT:C47824</a></li></ol></td><td><ol start="8"><li>Classification models</li><li>Naïve Bayesian model</li><li>AUC–ROC</li><li>5-fold cross validation</li><li>Smiles descriptors</li><li>malaria dataset</li></ol></td></tr><tr><td>Model</td><td><ol start="14"><li>bqmodel:isDescribedBy</li></ol><ol start="15"><li>bqmodel:isDescribedBy</li></ol><ol start="16"><li>bqmodel:isDescribedBy</li></ol><p><br></p></td><td><ol start="14"><li>Online Web server</li><li>Ersilia Model Hub</li><li>PubMed Identification Number PMID</li></ol></td><td><ol start="14"><li><a href="https://www.ebi.ac.uk/chembl/maip/https://github.com/ersilia-os/eos4zfyhttps://identifiers.org/pubmed:33618772">https://www.ebi.ac.uk/chembl/maip/<br></a></li><li><a href="https://www.ebi.ac.uk/chembl/maip/https://github.com/ersilia-os/eos4zfyhttps://identifiers.org/pubmed:33618772">https://github.com/ersilia-os/eos4zfy<br></a></li><li><a href="https://www.ebi.ac.uk/chembl/maip/https://github.com/ersilia-os/eos4zfyhttps://identifiers.org/pubmed:33618772">https://identifiers.org/pubmed:33618772</a></li></ol></td><td><ol start="14"><li>MAIP web platform - Source code</li><li>Ersilia Incorporation URL</li><li>PubMed URL</li></ol></td></tr></tbody></table>

### 7. Contextualize the Computational Metadata by adding DOME

The **DOME annotation** provides more contexts to the computational metadata by identifying which section of the modelling the metadata belong to.&#x20;

Adding DOME to the table shows this;&#x20;

| Metadata                                                                                                                                                                    | DOME                                                                                                                           |
| --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| <p>classification models</p><p>naïve Bayesian model</p><p>AUC–ROC</p><p>5-fold cross validation</p>                                                                         | <p>Optimization-Algorithm</p><p>Optimization-Algorithm</p><p>Evaluation-Performance Measure</p><p>Evaluation-Method</p>        |
| <p>MAIP web platform - Source code</p><p>Ersilia Incorporation URL</p><p>Smiles descriptors</p><p>malaria dataset</p><p>predictions of potential Antimalarial compounds</p> | <p>Model-Executable form</p><p>Model-Executable form</p><p>Data-Input</p><p>Data-Source</p><p>Model-Output; Classification</p> |

## &#x20;Resources & References

1. [The Use Case Published Annotation](https://www.ebi.ac.uk/biomodels/MODEL2405210002#Files)
2. [BioModels Annotation SOP](https://drive.google.com/file/d/1JqjcH0T0UTWMuBj-scIMwsyt2z38A0vp/view)
3. [Dome Annotation](https://www.nature.com/articles/s41592-021-01205-4)
4. [Ersilia Model for Use Case](https://github.com/ersilia-os/eos4zfy)
5. [Associated Publication for the Use Case](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-021-00487-2)
6. [Annotation Template](https://docs.google.com/spreadsheets/d/1bsCkN5Ugmo3tSF4oc-NGAQaZYfTP9LIABezrov-Vj-0/edit#gid=0)
7. [BioModelsML Publication](https://www.biorxiv.org/content/10.1101/2023.05.22.540599v1.full)
