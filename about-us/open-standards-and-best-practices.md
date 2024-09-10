---
description: This page describes Ersilia's Open Standards and Best Practices Principles
---

# Open Standards and Best Practices

## Open Standards

We follow [DPG Alliance tips for Open Standards](https://github.com/DPGAlliance/publicgoods-candidates/blob/main/help-center/open-standards.md). Below is a summary and proof of adherence to these standards.

| Concept                                                   | Comment                                                                                                                                                                                                                                                                      | Proof of Adherence                                                                                                                                                                                                                                                                                                                                                                                                                   |
| --------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Accessibility, Security, Authentication and Authorization | Ersilia uses GitHub primarily. Accessibility, security and authentication/authorization is therefore reliant on GitHub solutions. We use DependaBot to monitor package dependencies and identify security liabilities. GitHub Secrets are used intensively in our workflows. | See our [GitHub profile](https://github.com/ersilia-os) for more information. A summary of our GitHub security usage can be found [here](https://github.com/ersilia-os/ersilia/security).                                                                                                                                                                                                                                            |
| Internationalization                                      | UTF-8 encodings is used in our scripts, most of them written in Python. Code is formatted with Black.                                                                                                                                                                        | You can see our main codebase [here](https://github.com/ersilia-os/ersilia).                                                                                                                                                                                                                                                                                                                                                         |
| Application Programming Interfaces (APIs)                 | OpenAPI, especially via Swagger UI as facilitated by BentoML.                                                                                                                                                                                                                | Check one of our deployed models [here](https://eos80ch-m365k.ondigitalocean.app/).                                                                                                                                                                                                                                                                                                                                                  |
| Data Exchange and Configuration Formats                   | We primarily use YAML, JSON, CSV and TOML formats.                                                                                                                                                                                                                           | [Setup file](https://github.com/ersilia-os/ersilia/blob/master/pyproject.toml) in TOML format. [Metadata file](https://github.com/ersilia-os/eos3b5e/blob/main/metadata.json) in JSON. [Workflow file](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/pr\_check.yml) in YAML. [Data file](https://github.com/ersilia-os/pharmacogx-embeddings/blob/main/data/chemical\_descriptors/drug\_molecules.csv) in CSV. |
| Standard Content Formats and Multimedia                   | Content and multimedia are not our main assets. Internally, we store documents and media files with standard formats.                                                                                                                                                        | An example of a one-pager in PDF format can be found [here](https://drive.google.com/file/d/1Xxgpjh3gCQdD\_MqEDxweJIPY\_1JGKSIN/view?usp=sharing). We store videos in MP4 format, and upload them to [Youtube](https://www.youtube.com/channel/UCeioZf4Qj4hWi3O5Ta2k-xQ).                                                                                                                                                            |

## Best Practices and Principles

We follow [DPG Alliance tips for Best Practices and Principles](https://github.com/DPGAlliance/publicgoods-candidates/blob/main/help-center/best-practices.md). Below is a summary and proof of adherence to these principles.



| Concept                                                       | Comment                                                                                                                                   | Proof of Adherence                                                                                                                                                                                                                                                                                    |
| ------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| ICT4D                                                         | We endorse the [Principles for Digital Development](https://digitalprinciples.org).                                                       | Signed proof of endorsement available [here](https://github.com/ersilia-os/ersilia/security).                                                                                                                                                                                                         |
| User stories                                                  | We share community stories via our blog publication. We are active on social media and use GitHub Issues for the more technical aspects.  | Ersilia stories in [Medium](https://medium.com/ersiliaio). Activity in [GitHub Issues](https://github.com/ersilia-os/ersilia/issues). Other stories by the community: [GDI example](https://www.gooddatainstitute.com/post/pipeline-dreams-smiles-all-around-as-mlops-boosts-drug-discovery).         |
| Change management and version control                         | We use all version control functionalities of GitHub.                                                                                     | [Commit history](https://github.com/ersilia-os/ersilia/commits/master/) of our main repository.                                                                                                                                                                                                       |
| Test driven development using automated tests                 | Extensive and customized tests are performed, both for the main code and the models.                                                      | [Testing module](https://github.com/ersilia-os/ersilia/tree/master/test) in the Ersilia CLI. Test [workflow](https://github.com/ersilia-os/ersilia/actions/workflows/pr\_check.yml).                                                                                                                  |
| CI/CD                                                         | We use GitHub Actions for our CI/CD workflows, especially for model incorporation.                                                        | A subset of GitHub Actions can be found [here](https://github.com/ersilia-os/ersilia/actions/workflows/pr\_check.yml).                                                                                                                                                                                |
| Code review                                                   | We operate via pull requests and code reviewing is necessary previous to approval.                                                        | As an example, a video showing our model incorporation and reviewing process if found [here](https://youtu.be/I7dYI4ZF7Q0?si=3Ugo8zeqgBd5gQU-).                                                                                                                                                       |
| Agile development                                             | We use GitHub projects following the "Epic", "Sprint", etc. logic.                                                                        | [Agile template](https://github.com/ersilia-os/ersilia/blob/master/.github/ISSUE\_TEMPLATE/project-issue.yml) for defining projects.                                                                                                                                                                  |
| Modularity and Maintainability, Reusability and extensibility | The Ersilia Model Hub is highly modular, with each model corresponding to one artifact stored in an isolated environment.                 | [DockerHub](https://hub.docker.com/u/ersiliaos) model registry of Ersilia.                                                                                                                                                                                                                            |
| Component based architecture                                  | By design, our architecture is component based.                                                                                           | Find a high-level diagram [here](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/introduction).                                                                                                                                                                                             |
| Cloud Computing                                               | We have infrastructure as service scripts for AWS. However, the Ersilia Model Hub is designed to be cloud-agnostic.                       | Example of AWS integration [scripts](https://github.com/ersilia-os/model-inference-pipeline) for the Ersilia Model Hub.                                                                                                                                                                               |
| AI/ML                                                         | We use ONNX format for interoperability as much as possible in our assets.                                                                | This in-house library ([olinda](https://github.com/ersilia-os/olinda)) is specifically dedicated to distilling and converting models to ONNX.                                                                                                                                                         |
| Data Principles                                               | We follow FAIR data principles in all our projects. Most of our projects are collaborative.                                               | As an example, see our data sharing efforts as part of the [Ligand Discovery](https://github.com/ligand-discovery) project. These efforts were led by Ersilia.                                                                                                                                        |
| User Interface and User Experience (UI/UX)                    | Graphical user interfaces and a simple command-line interface are available.                                                              | [Ersilia CLI](https://github.com/ersilia-os/ersilia) is our primary tool. For some selected models, we offer online deployment with a [GUI](https://github.com/ersilia-os/ersilia-gui).                                                                                                               |
| Coding Styles and Standards                                   | We follow PEP 8 and format with Black.                                                                                                    | Black label in our [README file](https://github.com/ersilia-os/ersilia?tab=readme-ov-file).                                                                                                                                                                                                           |
| Open Source                                                   | We are a fully open source project and are part of several open source communities.                                                       | Software Sustainability Institute [cohort](https://www.software.ac.uk/blog/ersilia-ai-tool-drug-discovery-africa). Code for Science and Society [cohort](https://www.codeforsociety.org/eventfund/grantees/bringing-data-science-and-ai-ml-tools-to-infectious-disease-research-a-hands-on-workshop). |