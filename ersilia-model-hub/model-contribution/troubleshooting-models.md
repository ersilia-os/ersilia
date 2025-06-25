---
description: >-
  This page describes a few steps you can take when a Ersilia model is not
  working.
---

# Troubleshooting models

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
ersilia -v fetch eos9ei3 --from_dir eos9ei3 > out.log 2>&1
```

Always:

* Run the fetch command in verbose mode (`-v`) to print the output in the terminal
* Use the `--from_dir` flag to fetch from the locally cloned repository instead of its online version
* Copy the log to a file `out.log 2>&1`

We need to then carefully read the output that the fetch step is giving. We need to <mark style="color:red;">**read the whole fetch file**</mark>. In this log file, look for Error messages, package dependency incompatibilities, memory or connection issues, it can give you a good hint of what is happening. If it is an easy fix, you can go ahead an try it out by following the steps below.

If you encounter serious problems, please open a Bug Report issue in the Ersilia page and paste there as much information as possible, including the log file and which system you are using.

## 3. Test the model locally

#### 1. Create the conda environment

Since the fetch has failed, we do not have the necessary conda environment created. Depending on the step at which fetch is failing, you might have a conda environment with the model name, but it won't be complete. We recommend deleting it and starting anew, following the `install.yml` instructions.

In our example case, the `install.yml` looks like:

```yaml
python: "3.10"

commands:
  - ["pip", "rdkit", "2023.3.1"]
```

It seems the model runs on python 3.10 and it only requires the RDKIT package to work. If there is more than one package install required, please install them one by one to identify any possible package dependencies or deprecated packages.

```bash
conda create -n eos9ei3 python=3.10
conda activate eos9ei3
pip install rdkit==2023.3.1
```

#### 2. Run the model from run.sh

We have replicated the first steps of the Ersilia `fetch` command. If there are no package dependency issues and the conda environment is complete, we can go onto testing the model directly from the `run.sh` file. Remember to create a mock file for this purpose.

```bash
cd eos9ei3/model/framework
bash run.sh . examples/run_input.csv output.csv #compare with the actual run_output
```

This command will print an output on the terminal, which will give us hints of what can be the problem. Most typically, the errors are due to:

* Importing relative packages: if the import paths are not well specified, please reformat them to use absolute paths, avoiding future clashes.
* Input and output adapters: add print statements in the code to see that the input and output are in the right format, and modify them if needed.
* GPU - CPU issues: models that use pyTorch or other packages might have GPU specific configurations that need to be tweaked in order to work in most systems.

Once we are able to successfully run the run.sh model, we need to try it with Ersilia. Hopefully, the local fetch will now work. If it does not, go through the log file to obtain a hint of what might be going wrong or proceed to step 4 (advanced users).

```
ersilia -v fetch eos9ei3 --from_dir eos9ei3 > out.log 2>&1
```

## 4. Package the model locally

Our models are packaged using FastAPI via the [ersilia-pack](https://github.com/ersilia-os/ersilia-pack) implementation. If you are able to run the model from the run.sh file but unable to make it work within Ersilia, try packaging it locally to see where the source of the error is. We assume the model is cloned in the working directory and the conda environment of the model exists to be able to run the following code.

#### Install Ersilia Pack in its own environment

```bash
git clone https://github.com/ersilia-os/ersilia-pack
cd ersilia-pack
conda create -n ersiliapack python=3.12
conda activate ersiliapack
pip install -e .
cd ..
```

#### Bundle the model

Run the following command in your terminal, assuming you are in the root where the model directory is found

```bash
ersilia_model_pack --repo_path eos9ei3 --bundles_repo_path ~/eos/repository --conda_env_name eos9ei3 && ersilia_model_serve --bundle_path ~/eos/repository/eos9ei3 --port 45220
```

This will start the model package locally, with a message in the terminal like like:

```bash
INFO:     Will watch for changes in these directories: ['/home/ubuntu/eos/repository/eos4tcc/20250603-f8020c82-101e-4893-addf-8a85dc799323/app', '/home/ubuntu/eos/repository/eos4tcc/20250603-f8020c82-101e-4893-addf-8a85dc799323/model/framework']
INFO:     Uvicorn running on http://0.0.0.0:45220 (Press CTRL+C to quit)
INFO:     Started reloader process [67143] using WatchFiles
INFO:     Started server process [67145]
INFO:     Waiting for application startup.
Redis not connected
INFO:     Application startup complete.
```

Go to the URL specified by the INFO (http://0.0.0.0:45220). Do not close the terminal

#### Run the packaged model

The URL will bring you to the endpoints available for the model. You need to go to /docs to see them, so type in your browser: http://0.0.0.0:45220/docs

Go to the Run endpoint, select Try it out and run the example input provided. If it works, you'll see the expected output displayed. If it does not, you'll see an error message. Go back to the terminal and read through the exact error that happened when trying the Run endpoint (you need to scroll through the terminal to identify it)

#### Close the open port

Once you have completed the tests, is best to kill the process in the open port. For this you need to identify the Process ID (PID) and kill it:

```bash
lsof -i :45220
COMMAND   PID   USER   FD   TYPE DEVICE SIZE/OFF NODE NAME
python3  12345 username   12u  IPv4 123456      0t0  TCP *:45220 (LISTEN)
kill -9 12345
```

## 4. Update the model

* Remove any temporal edits, like `print` statements.
* Ensure the packages listed in the `install.yml` are updated to the version that is working.
* Revise the `.gitattributes` file.

When you are ready to push the changes to your fork of the model, open a PR to the main branch. As in the Model Incorporation, a series of [automated tests](example-of-the-model-incorporation-workflow.md#open-a-pull-request) will be triggered. Please check their result before moving on.

## 5. Check the Docker image

It may happen that workflows fail to build a Docker image, or that a Docker image is built and uploaded to DockerHub, but nonetheless the model does not successfully run. Below, we try to cover various scenarios and provide guidelines on how to troubleshoot them.

### Docker image is available in DockerHub, but it does not run successfully

If an image is already available in Ersilia's DockerHub but it fails to run successfully, then the best is to simply pull the image with the Docker CLI, and inspect it from the command line.

```bash
# docker pull ersiliaos/$MODEL_ID:[latest,dev-amd64,dev-arm64]
docker pull ersiliaos/eos5axz:latest
docker run -it --entrypoint /bin/bash ersiliaos/eos5axz:latest
```

Now your terminal should let you inspect the content of the docker container.

* If you model did not require any Conda installation, then your system Python should already have all dependencies installed. You can check it by simply typing `python` and trying to import any of the model packages (for example, `import rdkit`).
* If your model did require Conda installation, then in principle you should already be inside a default Conda environment with all your packages available.

To test the model inside the Docker container, the best is to go to the `framework` folder and simply run it using the `run.sh` file:

```bash
# cd /bundles/eos5axz/[random_id]/model/framework
cd /bundles/eos5axz/20250507-a54cd98a-fcb7-4989-9e99-1f10018ac004/model/framework
bash run.sh . examples/run_input.csv examples/run_output.csv
```

### Docker image works on AMD64, but it fails in ARM64

Oftentimes, building models on AMD64 platforms works but, unfortunately, they fail on ARM64. In that case, there can be several reasons for the failure:

* **Different packages between AMD64 and ARM64:** In that case, it will be necessary to specify the packages on a platform-specific way in the model repository. Broadly speaking, there are two scenarios:
  * Compilers are not available in ARM64
  * Model-specific packages have different versions/binaries in AMD64 vs ARM64
* **Different model code between AMD64 and ARM64:** This is a more rare scenario, but it can happen that code for running the models on ARM64 is different that that of AMD64. In that case, it will be necessary to provide code both for ARM64 and AMD64 separately.

Generally, the best way to troubleshoot ARM64-related issues is to build models from source in an ARM64 platform, including Apple M1/2/3 chips.

{% hint style="danger" %}
Note that we do not have a well-established way to specify platform-specific dependencies in the model repository, nor to specify platform-specific code. We are currently working on this.
{% endhint %}

{% hint style="warning" %}
ARM64 platforms (for example, Apple M1 chips) do not work with Python versions below 3.8. Please check that the Python version of your model is 3.8 or above. Otherwise, the model will not get packaged successfully for ARM64 architectures.
{% endhint %}

### Docker image was not pushed to DockerHub successfully

If your Docker image was not pushed to DockerHub successfully on the workflows, then you can try to install the model from source in an interactive way.

We have several base images containing the [Ersilia Pack](https://github.com/ersilia-os/ersilia-pack) repository. The naming convention is as follows:

* **ersiliaos/ersiliapack-py38** will contain Ersilia Pack and Python 3.8 as a system Python. Similarly, ersiliaos/ersiliapack-py312 will contain Python 3.12. You should consider this image if your model does _not_ have any Conda package.
* **ersiliaos/ersiliapack-conda-py38** will contain Ersilia Pack along with Conda (Python 3.8). You should consider this image if your model has Conda packages.

First, pull the image and run it in interactive mode:

```bash
docker pull ersiliaos/ersiliapack-py312:latest
docker run -it --entrypoint /bin/bash ersiliaos/ersiliapack-py312:latest
```

By default, you will be place in the  `/root`  directory in the container. You will notice that this directory is essentially empty.

You'll need to download the model source to work with it. The easiest is to get if from S3:

```bash
apt-get update
apt-get install wget
apt-get install unzip

wget https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos5axz.zip
unzip eos5axz.zip
```

Now you should have the model folder available. Ersilia is _not_ installed inside the base image, so you need to try and package the model with Ersilia Pack:

```bash
# ersilia_model_pack --repo_path $REPO_PATH --bundles_repo_path $BUNDLE_PATH
ersilia_model_pack --repo_path eos5axz --bundles_repo_path bundles
```

For more information, visit the [Ersilia Pack](https://github.com/ersilia-os/ersilia-pack) documentation.
