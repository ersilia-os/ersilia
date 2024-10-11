# Development
This file highlight the installation process of getting started with Ersilia, Ersilia Artifacts and Creating a PR.

## Getting Started
Ersilia can be run either on local computer or with the use of GitHub codespace. On your local computer, ersilia is designed to run specific environment
namely: **Linux and Linux-like environment** (MacOS, WSL).

- ### Installation of Ersilia on local computer
To get started with Ersilia installation, there are few dependencies that needs to be installed before proceeding to installation of Ersilia. Please check the package requirement
in the [Installation Guide](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/installation)
1. Once the required dependencies has been met, open a terminal and set up a Conda environment.

```
# create a conda environment
conda create -n ersilia python=3.10
# activate the environment
conda activate ersilia
```
2. Clone the repository and run to install the Ersilia Python package.
```
# clone from github
git clone https://github.com/ersilia-os/ersilia.git
cd ersilia
# install with pip (use -e for developer mode)
pip install -e .
```
3. To test that the CLI works, explore this command
```
# output to see command options for ersilia
ersilia --help
```
4. Once ersilia is recognized in the CLI, test any model of your choice.
- First, you can check the list of available model in ersilia catalog
```
#see ersilia's model catalog
ersilia catalog
```
- Next, you can run the model
```
ersilia -v fetch eos2ta5 #fetch a model
ersilia serve eos2ta5 #serve a model
ersilia -v api run -i "CCCC" #run a model
```

Note: The usage -v (verbose) flag is of utmost importance for developers as it provides detailed logs which helps to debug and troubleshoot issues.

- ### Using GitHub Codespace
If difficulties are encountered when setting up Ersilia locally, another cloud-based development environment which suppoort linux is a great alternative.
GitHub codespaces allows you to run and develop Ersilia directly from your browser.

To use Github Codespaces:

 **1. Fork and Clone the Repository**
 
 - Navigate [Ersilia repository](https://github.com/ersilia-os/ersilia) on GitHub.
 - Click on the **Fork** button in the top right corner to create a fork.
 - Go to the forked repository, Click on the <> **Code**.
 - Click **Create codespace**.
   
A pop up page "Setting up your codespace" is displayed on the screen. After a while, the codespace is created in a browser-based version using the built-in 
Visual studio code interface 

**2. Setting Up Codespace Environment**

- Open the terminal in codespace. To ensure you are in the root directory
`cd ersilia`
- Run this command
`pip install -e .`

**3. Run Ersilia Model**
- Directly from the code space terminal, you can fetch, serve and run Ersilia model.
To run a model with verbose output:
`ersilia -v run eos2ta5`


## Ersilia Artifacts
After successful installation of ersilia, the system creates an EOS directory in your `$HOME` path. This directory stores
ersilia artifacts that the local computer uses when doing any operation (model execution) in Ersilia.
Folders found in the `$HOME` directory includes:

**Folders**
- dest/: stores output generated when running a model.
- isaura/: cache model prediction which stores them in the local computer
- repository/: main directory to cache models files and datasets. This prevent redundant download by storing previously fetched model
- sessions/: Tracks user session
- tmp/: stores temporary file during model execution.


  



