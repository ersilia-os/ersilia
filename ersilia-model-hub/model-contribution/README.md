---
description: Documentation for those who want to contribute new models to the Hub
---

# Model contribution

## Setting up

To contribute models at Ersilia you need to be familiar with:

* GitHub
* Conda
* Python
* Docker

### Installation

It is recommended to install Ersilia locally in \[test] mode to enable for model testing upon incorporation.&#x20;

```bash
conda create -n ersilia python=3.12
conda activate ersilia
pip install ersilia[test]
```

Ersilia is only maintained for Linux/MacOS systems. If you are working with Windows, please use a Windows Substsystem Linux.

#### Additional installations

* DockerHub: all models incorporated in Ersilia will be Dockerized for easy deployment. While the dockerization step happens in the cloud (GitHub Actions) it is recommended to install docker for model testing purposes.
* Git-LFS: many model’s checkpoints are too large (>100MB) for GitHub storage. Make sure that [git-lfs](https://git-lfs.com/) is installed and active in your system to push large files to the model repository.&#x20;

## TL:DR

Summary of the contributor’s pipeline:

1. Read about and understand your model of choice if you are not its main developer.
2. Follow the original author’s instructions to run the model locally and note down the steps and requirements needed.
3. Open a Model Request issue under the Ersilia repository.
4. Fork the repository created upon approval of the Model Request by Ersilia maintainers.
5. Clone the model repository in your local and focus on the following:
   1. Include the model checkpoints in the `/checkpoints` directory. Track them as git-lfs files using the `.gitattributes` file if necessary
   2. Include any code needed to run the model as .py script files in the `/framework/code` directory
   3. Adapt the functions in `main.py` to actually run your model and provide the output in the specified format
   4. Create the `run_columns.csv` file including information about the model output.&#x20;
   5. Include the dependencies (specifying the version) in the `install.yml` file
   6. Complete all fields on the `metadata.yml`. Note that many are pre-defined, ensure you are using the correct notation.
6. Open a PR to the main branch of the model repository. A series of automated tests will be triggered. Ensure they are passing.
7. Check if your PR has been successfully merged and the model is uploaded to AWS S3 and DockerHub.
8. Delete your fork of the model repository.

If it is your first time incorporating a model in the Ersilia Model Hub, make sure to read with detail the following sections:

* [Model template](model-template.md): detailed explanation of the folder structure and all the files present in the eos-template.
* [Model incorporation workflow](example-of-the-model-incorporation-workflow.md): exemplified step by step model incorporation.
* [Troubleshooting models](troubleshooting-models.md): pointers of what to check if your incorporated model is not working.
