---
description: >-
  This page outlines our goal to create the reference database of purchasable
  compound collections with potential interest to the field of AMR
---

# AMR chemical collections

## Background

We have plans to develop the **Ersilia AMR Chemical Collections (E-AMR-CC)** web app, initially focused on _Klebsiella pneumoniae_. E-AMR-CC will consist of a catalog of bespoke chemical libraries with predicted anti-_Klebsiella_ activity. The resource will be extensible to other species in the future.

## Current repositories

Below is a table of the current code repositories available at Ersilia that will be key to develop the E-AMR-CC resource.

| Repository Name                                                                                          | Description                                                                                                            | Relevance to E-AMR-CC                                                                                                                                           |
| -------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [Ersilia Model Hub](https://github.com/ersilia-os/ersilia)                                               | Main CLI tool from Ersilia, offering easy-to-use functionality to fetch and serve AI/ML models for drug discovery.     | The Ersilia Model Hub contains general-purpose models for drug discovery, which will be used to annotate the compound libraries.                                |
| [ZairaChem](https://github.com/ersilia-os/zaira-chem)                                                    | End-to-end modeling based on ensembles of descriptors and autoML tools.                                                | For the supervised AI/ML models (phenotypic), this will be the main pipeline to build and distill models.                                                       |
| [ChEMBL Tasks for Fine Tuning](https://github.com/ersilia-os/chembl-tasks-for-finetuning)                | Generic data collection tool from ChEMBL, with a focus on AI/ML-readiness.                                             | This tool will be used to fine-tune TabPFN models, used especially in the modellability assessment.                                                             |
| [ChEMBL Binary Tasks](https://github.com/ersilia-os/chembl-binary-tasks)                                 | Workflow to collect antimicrobial data from ChEMBL, including data preparation with LLMs, and modellability assesment. | The workflow has been tested on _A. baumannii_. It will now be used on _Klebsiella pneumoniae_.                                                                 |
| [_A. baumannii_ Enamine REAL screening](https://github.com/ersilia-os/abaumannii-enamine-real-screening) | Screening of the Enamine REAL 1B Lead-like library against 30 _A. baumannii_ models.                                   | This is an example repository demonstrating that billion-scale screening is feasible for anti-_Klebsiella_ activity prediction models.                          |
| [_Mtb_ Targeted Protein Degradation](https://github.com/ersilia-os/mtb-targeted-protein-degradation)     | GC-ADDA4TB project aimed at discovering pan-engaging warheads against essential tRNA synthetases.                      | Structure-based approach containing robust detection of binding pockets in ensembles of structures. A similar strategy will be applied to _Klebsiella_ targets. |

