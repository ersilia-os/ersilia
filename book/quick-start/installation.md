# Installation

## Before you start

There are a few third-party **pre-requisites** and you need to have them installed on your computer before proceeding with the installation of the Ersilia tool. Don't worry, these dependencies are very popular, well documented and continuously maintained. Most probably you know them already:

### Python and Conda

We have tested our tools on **Python 3.7** and above. Please make sure you have the right Python installation on your computer. Visit the [official Python site](https://www.python.org) to learn more.

Although _not_ strictly necessary, we encourage the use of **Conda environments**. If Conda is available in the local device, the Ersilia tool will try to use it in order to handle dependencies safely and efficiently. Conda can be installed following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Many people choose to use [Anaconda](https://docs.anaconda.com/anaconda/install/index.html) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installers.

### Docker

Docker containers are an excellent way to share applications and simplify the management of system requirements and configurations. Please [install Docker](https://www.docker.com) to ensure that all our AI/ML assets will work smoothly in your local device.

### Git and Git LFS

All of our code is hosted in the [Ersilia Github Organization](https://github.com/ersilia-os) profile. The Ersilia tool needs to fetch code from GitHub, so make sure Git is installed on your device. You can follow [these instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git). In addition, we highly recommend the fantastic [GitHub CLI](https://cli.github.com/manual/installation), which is used by Ersilia if available.

Some pre-trained models, especially the ones that contain many parameters, require **Git Large File Storage** (LFS). Make sure this Git extension is [installed](https://git-lfs.github.com) before proceeding.

## Install on Linux and MacOSX

Open a terminal. The best is to set up a [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) environment.

```bash
# create a conda environment
conda create -n ersilia python=3.7
# activate the environment
conda activate ersilia
```

Then, simply install the Ersilia Python package.

```bash
# clone from github
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
# install with pip (use -e for developer mode)
pip install -e .
```

You should be done! Quickly check that the CLI works on your terminal.

```bash
# see ersilia CLI options
ersilia --help
```

### The Isaura data lake

We highly recommend installation of the [Isaura](https://github.com/ersilia-os/isaura) data lake. With Isaura, you will be able to cache your model predictions (i.e. store them in your local computer). Isaura is a relatively light Python package:

```
# clone from github
git clone https://github.com/ersilia-os/isaura.git
cd isaura
# install with pip (use -e for developer mode)
pip install -e .
```

## Install on Windows

{% hint style="warning" %}
We are not testing Windows installation consistently. If you encounter problems, please reach out to us at hello@ersilia.io.
{% endhint %}

We recommend that you install Ubuntu on your Windows machine. This can be now done very easily with WSL. You will need at least Windows 10.

Open a Power Shell with Admin permissions and type:

```
wsl --install
```

Then simply [install the Ubuntu terminal](https://www.microsoft.com/en-us/p/ubuntu/9nblggh4msv6#activetab=pivot:overviewtab) on Windows.

Inside the Ubuntu terminal, you can now follow the installation instructions [above](installation.md#install-on-linux-and-macosx).

## Command-line only installation snippets for Ubuntu

Below you can find a few snippets that can help you install dependencies in a Ubuntu terminal.

### The gcc compiler

Probably you have the gcc compiler installed already. Just in case:

```
sudo apt install build-essential
```

### Conda package manager

If you don't have the Conda package manager yet, we suggest you install Miniconda:

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
~/miniconda3/bin/conda init zsh
```

### GitHub CLI

Once Conda is installed (see above), you can use it to install the fantastic GitHub CLI:

```
conda install gh -c conda-forge
```

Use the GitHub CLI to login. This may be helpful if you have contributor permissions at Ersilia. Type:

```
gh auth login
```

And then follow the instructions.

### Git LFS

Git Large File Storage (LFS) can be installed from Conda as well:

```
conda install git-lfs -c conda-forge
```

Activate Git LFS:

```
git-lfs install
```
