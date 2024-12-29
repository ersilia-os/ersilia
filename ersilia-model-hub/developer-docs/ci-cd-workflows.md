---
description: >-
  Ersilia relies heavily on GitHub Actions for CI/CD. Here is a high-level
  summary of workflows involved in maintaining the Ersilia CLI and code-bases
  from all  the Ersilia maintained models.
---

# CI/CD workflows

## Ersilia repository

These workflows ensure that we push quality code to Ersilia, incorporating testing by individual components as well as end-to-end.&#x20;

### Testing

#### Testing Ersilia on a Pull Request

* Objective: ensure Ersilia continues to work with any new code merging from contributors and maintainers to the master branch
* Workflow name: Deploy and test ersilia on PR
* Workflow file: [pr\_check.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/pr_check.yml)
* Jobs:&#x20;
  * `Build`: creates a conda environment (py3.10 currently) and installs Ersilia from the open PR. Then it tests the following functions: catalog and model running (fetch, serve, info and run) from GitHub, S3 and DockerHub. A single molecule is passed as input, and the output is analysed <mark style="background-color:red;">...</mark>
* Run: automatically when a PR is open on the Ersilia master branch

#### Testing and cleaning up

* Objective:
* Workflow name: Ersilia tests, installation checks, and cleanup of model request template
* Workflow file: [tests\_and\_cleanup.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/tests_and_cleanup.yml)
* Jobs:
  * `install-ersilia`, installs Ersilia from Python 3.8 to Python 3.12 on an Ubuntu base.
  * `test-docker`  is the next step to follow after we confirm that Ersilia can be installed in environments created with different Python versions. This step tests whether dockerized models can be fetched using the ersilia CLI, ie the happy flow works without issues.
  * `test-os`  Since Ersilia is supported only on Linux distributions, and MacOS, this job ensures that Ersilia can be installed and run on both of these platforms, presently only utilizing a Python 3.10 environment.&#x20;
  * `run-pytest`  This is the main step that carries out unit testing within Ersilia ensuring that different components continue to work without issues.
  * `run-cli-test-single` This an integration level test that ensures different components within Ersilia interact cohesively in running a single model using the Ersilia CLI. More details on how this is implemented can be found [here](testing-playground.md).
  * `run-cli-test-multiple` Similar to the job described above, this job tests the functionality in Ersilia to simultaneously run multiple models. More details around implementation of the specific tests and the testing frameworks used can be found [here](testing-playground.md).
  * `test-colab-notebook` This test ensures that Ersilia can run in IPython notebooks, such as Jupyter or Google Colab. Specifically, this test checks whether a notebook containing Ersilia code can be executed successfully.
  * `update-model-request-template` This job within the workflow only runs when the PR or the commit on `master` branch change the file containing [model tags](https://github.com/ersilia-os/ersilia/blob/master/ersilia/hub/content/metadata/tag.txt). We expect model tags to be completely sorted because they are used to populate the [Model Request Issue](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&projects=\&template=model_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) (more on this below). The job sorts the tags file if it is changed, and then updates the Model Request Issue [template](https://github.com/ersilia-os/ersilia/blob/master/.github/ISSUE_TEMPLATE/model_request.yml) with the changes, and commits them back to the master branch.
* Run:

### Packaging

#### Creating Docker image of Ersilia's CLI

* Objective: package Ersilia into an image to use as a base image to containerise individual models as well as in development containers.
* Workflow name: Upload Ersilia base image to DockerHub
* Workflow file: [ersilia-base-image-to-dockerhub.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/ersilia-base-image-to-dockerhub.yml)
* Jobs:
  * `upload_ersilia_base_to_dockerhub` builds two versions of the Ersilia image,  a lighter version tagged as _`ersiliaos/base-multistage-condapack` ,_ and a legacy much bloated version tagged as _`ersiliaos/base-legacy-bentoml` ._&#x20;
* Run: configured to trigger in multiple ways - whenever code is pushed on to the `master` branch, when a new release is created from Ersilia at the beginning of each month (more on this below), and also manually by developers and maintainers.&#x20;

#### Ersilia Release

* Objective: Ersilia is distributed as source, as a Docker image, as a PyPI package, and as a Conda Package.&#x20;
* Workflow name: Ersilia Release
* Workflow file: [publish.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/publish.yaml)
* Jobs:
  * `version` This is the job that increments the version of the code-base, and creates a new release tag, while also updating this information in the [CITATION.cff](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff), and [CODEMETA.json](https://github.com/ersilia-os/ersilia/blob/master/codemeta.json) files.&#x20;
  * `gh-release` Utilizing the tag created in the previous job, a GitHub release is created in this job of the workflow. The GitHub CLI is used to automatically create this release and utilize the commit history from the `master` branch to generate release notes.
  * `pypi-release` Finally, we utilize the release artifacts created from the GitHub release job of the workflow to publish a new version of Ersilia on the Python Packaging Index (PyPI).
  * The completion of this release workflow triggers the Docker build workflow mentioned above.
  *   `Conda-Forge Release`

      Ersilia is released on conda-forge using by way of updating its [feedstock](https://github.com/conda-forge/ersilia-feedstock) maintained within the conda-forge organization on GitHub. Upon the creation of a new GitHub release, an automated pull request is created in this feedstock repository which is then manually reviewed and merged by an Ersilia maintainer.
* Run: At the beginning of each month, a release [workflow](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/publish.yaml) is executed at 3:00 AM UTC publishing a new release for Ersilia using the code present on the `master` branch.&#x20;

### Community

#### New Model Request

* Objective
* Workflow name:
* Workflow file:
* Jobs:
* Run:

## Model incorporation

### Testing

* XXX
* XXX

### Packaging

* XXX
* XXX

### Community

* XXX
* XXX

## Production and serving

* Ersilia Statistics
* Ersilia Maintenance
