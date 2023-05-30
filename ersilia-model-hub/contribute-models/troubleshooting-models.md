---
description: >-
  This page describes a few steps you can take when a Ersilia model is not
  working.
---

# Troubleshooting models

{% hint style="info" %}
This information is orientative. Each model will have its own specific requirements. If it is your first time contributing to the Ersilia Model Hub, please get acquainted with its structure before tackling these issues.
{% endhint %}

{% hint style="warning" %}
This documentation refers only to technical issues, package dependencies and the similar. We are not focusing on whether a model predicts more or less accurate predictions for a specific use case.
{% endhint %}

When Ersilia fetches a model from our online repository it automatically tests it using a random 3-molecule input. If the model is unable to produce an output, you will receive the following error and the model will not be fetched:

```python
"Error occured while fetching model: eos0abc"
```

If that is the case, go through the following steps:

## 1. Sanity checks

Make sure to:

1. Use the **latest version** of Ersilia. In your local clone of the Ersilia repository, please run `git pull --rebase`. If you installed Ersilia in the ersilia conda environment using `pip install -e .` you don't need to reinstall it, conda will use the updated codebase. If you installed it by simply running a `pip install` command, please reinstall Ersilia in your conda environment.
2. Install and activate `git-lfs`. The **Large File Storage** system from GitHub allows to store files > 500 MB in GitHub. It is possible that `git-lfs` is installed in your system but not activated. Typically, Ersilia will download the model from our AWS S3 storage, but some essential files might be stored in `git-lfs`.
3. Check that the `eos` folder is being created in your Home directory. There, a `dest` folder with the model subfolders should appear, this is the directory where Ersilia will clone the model. If the `eos` folder does not exist, you might need to revise your **User permissions for creating new folders**.
4. Identify if there is a **connection error.** Run the fetch command in verbose mode: `ersilia -v fetch eos0abc` and scan the output printed in the terminal for errors such as `TimeoutError: [Errno 60] Operation timed out` or `NewConnectionError`. Typically you will receive errors from the `urllib3` library if there is an internet issue. Please change network and make sure you are not using a VPN or firewall before trying again.
5. Make sure you have **enough free space** in your system. If running the fetch command in verbose model prints in the output a message like `No Space Left on Device`, free it up before proceeding. We recommend at least 10 GB of free disk space.

Once you have completed the above steps, try fetching the model one more time. If the problem persists, <mark style="color:red;">**do not insist**</mark> on model fetching from the online repository, follow the instructions below. We provide a step by step suggestion of what actions a contributor can take to identify the source of the problem. We are using the model from the [Example Workflow](example-of-the-model-incorporation-workflow.md) to run this debugging demo.

## 2. Fetch the model from local

The first step is to download the model in your system and use the fetch from local path option. This will avoid connection issues and failures in downloading from git-lfs or S3, as well as speed up the debugging process. Fork the repository and clone it.

```bash
git clone https://github.com/contributor-user/eos9ei3.git
ersilia -v fetch eos9ei3 --repo_path eos9ei3 > out.log 2>&1
```

Always:

* Run the fetch command in verbose mode (`-v`) to print the output in the terminal
* Use the `--repo_path` flag to fetch from the locally cloned repository instead of its online version
* Copy the log to a file `out.log 2>&1`

We need to then carefully read the output that the fetch step is giving. At the end, it typically says:&#x20;

```
Model API eos9ei3:predict did not produce an output
```

This is not the final error, it is simply stating that it could not calculate the test molecules, hence the empty output. We need to <mark style="color:red;">**read the whole fetch file**</mark>. In this log file, look for Error messages, package dependencies, memory or connection issues, it can give you a good hint of what is happening. If it is an easy fix, you can go ahead an try it out by following the steps below.

If you encounter serious problems, please open a Bug Report issue in the Ersilia page and paste there as much information as possible, including the log file and which system you are using.

## 3. Test the model locally

#### 1. Create the conda environment

Since the fetch has failed, we do not have the necessary conda environment created. Depending on the step at which fetch is failing, you might have a conda environment with the model name, but it won't be complete. We recommend deleting it and starting anew, following the Dockerfile instructions.

In our example case, the Dockerfile looks like:

```docker
FROM bentoml/model-server:0.11.0-py37
MAINTAINER ersilia

RUN conda install -c conda-forge rdkit=2021.03.4

WORKDIR /repo
COPY . /repo
```

It seems the model runs on python 3.7 and it only requires the RDKIT package to work. If there is more than one package install required, please install them one by one to identify any possible package dependencies or deprecated packages.

```bash
conda create -n eos9ei3 python=3.7
conda activate eos9ei3
pip install rdkit
```

#### 2. Run the model from run.sh

We have replicated the first steps of the Ersilia `fetch` command. If there are no package dependency issues and the conda environment is complete, we can go onto testing the model directly from the `run.sh` file. Remember to create a mock file for this purpose.

```bash
cd eos9ei3/model/framework
bash run.sh . ~/Desktop/test.csv ~/Desktop/output.csv
```

This command will print an output on the terminal, which will give us hints of what can be the problem. Most typically, the errors are due to:

* Importing relative packages: if the import paths are not well specified, please reformat them to use absolute paths, avoiding future clashes.
* Input and output adapters: add print statements in the code to see that the input and output are in the right format, and modify them if needed.
* GPU - CPU issues: models that use pyTorch or other packages might have GPU specific configurations that need to be tweaked in order to work in most systems.

Once we are able to successfully run the run.sh model, we need to try it with Ersilia. Hopefully, the local fetch will now work. If it does not, go through the log file to obtain a hint of what might be going wrong.

```
ersilia -v fetch eos9ei3 --repo_path eos9ei3 > out.log 2>&1
```

## 3. Update the model

* Remove any temporal edits, like `print` statements
* Ensure the packages listed in the `Dockerfile` are updated to the version that is working
* Revise the `.gitattributes` file

You are ready to push the changes to your fork of the model, and then open a PR to the main branch. As in the Model Incorporation, a series of [automated tests](example-of-the-model-incorporation-workflow.md#open-a-pull-request) will be triggered. Please check their result before moving on.
