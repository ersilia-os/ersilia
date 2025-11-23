# Development
This file highlight the installation process of getting started with Ersilia, Ersilia Artifacts and Creating a Pull Request.

## Getting Started
Ersilia can be run either on local computer or with the use of GitHub codespaces. Ersilia is designed to run on **Linux or Linux-like environments** namely: MacOS, Windows Sub-system for Linux (WSL) and at present is not supported on Windows systems.

- ### Installation of Ersilia on local computer
To get started with Ersilia installation, there are few dependencies that needs to be installed before proceeding to installation of Ersilia. Please, check the package requirement in the [Installation Guide](https://ersilia.gitbook.io/ersilia-book/ersilia-model-hub/installation)
1. Once the required dependencies has been installed, open a terminal and set up a Conda environment.

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
pip install -e .[test]
```
**The -e (editable) for developer mode ensures that whenever a user change the code, they don't need to reinstall ersilia for those changes to take effect.**

3. To test that the CLI works, explore this command
```
# output to see command options for ersilia
ersilia --help
```
4. Once ersilia is recognized in the CLI, test any model of your choice.
- First, you can check the list of available model in ersilia's catalog
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

- ### Using GitHub Codespaces
If difficulties are encountered while setting up Ersilia locally, another cloud-based development environment which suppoort linux is a great alternative.
GitHub codespaces allows you to run and develop Ersilia directly from your browser.

To use Github Codespaces:

 **1. Fork and Clone the Repository**
 
 - Navigate [Ersilia repository](https://github.com/ersilia-os/ersilia) on GitHub.
 - Click on the **Fork** button in the top right corner to create a fork.
 - Go to the forked repository, Click on the <> **Code**.
 - Click **Create codespace**.
   
A pop up page "Setting up your codespace" is displayed on the screen. After a while, the codespace is created in a browser-based version using the built-in Visual studio code interface 

**2. Setting Up GitHub Codespaces Environment**

- Once you open the codespaces environment, ersilia has been installed and ready for usage. To double check, you can verify by running this command `ersilia --version`. This
will display the version of ersilia installed.

**3. Run Ersilia Model**
- Directly from the code space terminal, you can fetch, serve and run Ersilia model.
To run a model with verbose output:
`ersilia -v run eos2ta5`

## Ersilia Artifacts
After successful installation of ersilia in the user's local computer, the system creates an EOS directory in the user's `$HOME` path. This directory stores ersilia artifacts that the local computer uses when doing any operation (model execution) in Ersilia.
Folders found in the `$HOME` directory includes:

**Folders**
- dest/: stores output generated when a model.
- repository/: main directory to cache models files and datasets. This prevent redundant download by storing previously fetched model
- sessions/: Tracks user session
- tmp/: stores temporary file during model execution.

## Creating a Pull Request
To get started, please read and follow our [Code of Conduct](https://github.com/ersilia-os/ersilia/blob/master/CODE_OF_CONDUCT.md).

Before submitting a pull request, please ensure your code has been tested. After making any changes, **run the test suites** to test the code properly. 

To do this, use this command to run the test with -v (verbose) output: `pytest -v`

**1. Via GitHub Codespaces Environment**
- In the codespace terminal, create a new branch
```
git checkout -b your-branch-name
```
- Make all necessary changes and commit them
```
git add .
git commit -m "Description of changes made"
```
- Push changes to GitHub
```
git push origin your-branch-name
```
- Click on **Create pull request** in the prompt that will be shown

**2. Via GitHub Desktop**
- Fork the repository
- Create a new branch and name it.
- Commit your changes
- Navigate the repository,Click on **Compare & pull request** 
- Set the main branch to **master** (name of ersilia repository main branch) and set the compare branch to your-branch-name
- Fill out the title section and give a brief description of any changes made
- Click **Create Pull Request**

  



