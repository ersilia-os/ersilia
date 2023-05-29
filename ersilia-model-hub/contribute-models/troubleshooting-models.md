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

1. Use the latest version of Ersilia. In your local clone of the Ersilia repository, please run `git pull --rebase`. If you installed Ersilia in the ersilia conda environment using `pip install -e .` you don't need to reinstall it, conda will use the updated codebase. If you installed it by simply running a `pip install` command, please reinstall Ersilia in your conda environment.
2. Install and activate `git-lfs`. The Large File Storage system from GitHub allows to store files > 500 MB in GitHub. It is possible that `git-lfs` is installed in your system but not activated. Typically, Ersilia will download the model from our AWS S3 storage, but some essential files might be stored in `git-lfs`.
3. Check that the `eos` folder is being created in your Home directory. There, a `dest` folder with the model subfolders should appear, this is the directory where Ersilia will clone the model. If the `eos` folder does not exist, you might need to revise your User permissions for creating new folders.
4. Identify if there is a connection error. Run the fetch command in verbose mode: `ersilia -v fetch eos0abc` and scan the output printed in the terminal for errors such as `TimeoutError: [Errno 60] Operation timed out` or `NewConnectionError`. Typically you will receive errors from the `urllib3` library if there is an internet issue. Please change network and make sure you are not using a VPN or firewall before trying again.
5. Make sure you have enough free space in your system. If running the fetch command in verbose model prints in the output a message like `No Space Left on Device`, free it up before proceeding. We recommend at least 10 GB of free disk space.

Once you have completed the above steps, try fetching the model one more time. If the problem persists, <mark style="color:red;">**do not insist**</mark> on model fetching from the online repository, follow the instructions below

## 2. Look at the fetch information

The first thing we need to do is carefully read the output that the fetch step is giving. At the end, it typically says&#x20;

```
Model API eos1af5:predict did not produce an output
```

This is not the final error, it is simply stating that it could not calculate the test molecules, hence the empty output. <mark style="color:red;">**We need to read the whole fetch file**</mark>. You can print it on the terminal using the verbose (-v) flag, or, even better, by copying the output to a log file:

```bash
ersilia -v fetch eos9ei3 >> out.log 2>&1
```

In this log file, look for Error messages, package dependencies, memory or connection issues, it can give you a good hint of what is happening. If it is an easy fix, you can go ahead an try it out by following the steps below.

If you encounter serious problems, please open a Bug Report issue in the Ersilia page and paste there as much information as possible, including the log file and which system you are using.

## 3. Test the model locally

Here, we provide a step by step suggestion of what actions a contributor can take to identify the source of the problem. We are using the model from the [Example Workflow](example-of-the-model-incorporation-workflow.md) to run this debugging demo.

### 1. Clone the model repository

The first thing we will do is clone the model in our local system, so we can avoid connection issues and failures in downloading from git-lfs or S3 as well as speeding up the debugging process.&#x20;

```bash
git clone https://github.com/ersilia-os/eos9ei3.git
cd eos9ei3
```

### 2. Create a conda environment

Since we have locally cloned the repository, we do not have the necessary conda environment to run it. We can manually create it by using the Dockerfile instructions. We recommend running each install line one by one to identify any possible package dependencies or deprecated packages.

```bash
conda create -n eos9ei3 python=3.7
conda activate eos9ei3
pip install rdkit
```

In this case, we only have one package specified, so we go ahead and pip install it.&#x20;

### 3. Run the model from run.sh

We now can try to run the model directly from the run.sh file. Remember to create a mock file for this purpose.

```bash
cd eos9ei3/model/framework
bash run.sh . ~/Desktop/test.csv ~/Desktop/output.csv
```

This will probably give us the information we need to debug the error.&#x20;

## 3. Update the model

