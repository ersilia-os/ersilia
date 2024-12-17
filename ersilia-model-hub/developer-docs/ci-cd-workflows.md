---
description: A summary of the
---

# CI/CD workflows

Ersilia relies heavily on GitHub Actions for CI/CD. Here is a high-level summary of workflows involved in&#x20;

## The main Ersilia repository

These workflows ensure that we push quality code to Ersilia, incorporating testing both at component level and tests for the end-to-end happy flows. Ersilia presently supports Python versions from 3.8 to 3.12, and only runs on Linux distributions and Mac OS. Read on to find out how we ensure these checks are carried out automatically.

### Testing

*   **Deploying and testing Ersilia on a Pull Request.**

    This [workflow](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/pr_check.yml) ensures that happy flows within Ersilia continue to work on every pull request that is created with the intention to push new code into Ersilia's master branch. The workflow, containing a single job called `build`, imitates steps that any Ersilia user would carry out - ie, it creates a new conda environment, currently only tied to Python 3.10, and installs ersilia from the user's PR in this environment. Once that is setup, the ersilia CLI is tested in multiple ways, such as displaying the catalog of Ersilia maintained models within the Ersilia Model Hub, as well as fetching, serving, and running a simple model from multiple sources.&#x20;
*   **Testing and Clean up**

    We have another [layer](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/tests_and_cleanup.yml) of testing that incorporates unit, and integration testing of the various components that make up the Ersilia code base. This workflow has multiple jobs, which are explained here:&#x20;

    * `install-ersilia`, which checks the feasibility of Ersilia installation within versions of Python 3.8 to Python 3.12 on an Ubuntu base.
    * `test-docker`  is the next step to follow after we confirm that Ersilia can be installed in environments created with different Python versions. This step tests whether dockerized models can be fetched using the ersilia CLI, ie the happy flow works without issues.
    * `test-os`  Since Ersilia is supported only on Linux distributions, and MacOS, this job ensures that Ersilia can be installed and run on both of these platforms, presently only utilizing a Python 3.10 environment.&#x20;
    * `run-pytest`  This is the main step that carries out unit testing within Ersilia ensuring that different components continue to work without issues.
    * `run-cli-test-single` This an integration level test that ensures different components within Ersilia interact cohesively in running a single model using the Ersilia CLI. More details on how this is implemented can be found [here](testing-playground.md).
    * `run-cli-test-multiple` Similar to the job described above, this job tests the functionality in Ersilia to simultaneously run multiple models. More details around implementation of the specific tests and the testing frameworks used can be found [here](testing-playground.md).
    * `test-colab-notebook` This test ensures that Ersilia can run in IPython notebooks, such as Jupyter or Google Colab. Specifically, this test checks whether a notebook containing Ersilia code can be executed successfully.
    * `update-model-request-template` This job within the workflow only runs when the PR or the commit on `master` branch change the file containing [model tags](https://github.com/ersilia-os/ersilia/blob/master/ersilia/hub/content/metadata/tag.txt). We expect model tags to be completely sorted because they are used to populate the [Model Request Issue](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&projects=\&template=model_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) (more on this below). The job sorts the tags file if it is changed, and then updates the Model Request Issue [template](https://github.com/ersilia-os/ersilia/blob/master/.github/ISSUE_TEMPLATE/model_request.yml) with the changes, and commits them back to the master branch.

### Packaging

*   **Creating Docker Image of the Ersilia CLI.**

    We heavily make use of Docker images at Ersilia, from both packaging Ersilia into an image to use in Development containers, to also using it as a base image to further containerize models. This [workflow](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/ersilia-base-image-to-dockerhub.yml) containing a single job called `upload_ersilia_base_to_dockerhub` is configured to trigger in multiple ways - whenever code is pushed on to the `master` branch, when a new release is created from Ersilia at the beginning of each month (more on this below), and also manually by developers and maintainers.&#x20;

    This build process supports building two versions of the image - a lighter version tagged as _`ersiliaos/base-multistage-condapack` ,_ and a legacy much bloated version tagged as _`ersiliaos/base-legacy-bentoml` ._&#x20;
*   **Ersilia Release**

    Ersilia is distributed as source, as a Docker image, as a PyPI package, and as a Conda Package. At the beginning of each month, a release [workflow](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/publish.yaml) is executed at 3:00 AM UTC publishing a new release for Ersilia using the code present on the `master` branch. This workflow has three main jobs:

    * `version` This is the job that increments the version of the code-base, and creates a new release tag, while also updating this information in the [CITATION.cff](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff), and [CODEMETA.json](https://github.com/ersilia-os/ersilia/blob/master/codemeta.json) files.&#x20;
    * `gh-release` Utilizing the tag created in the previous job, a GitHub release is created in this job of the workflow. The GitHub CLI is used to automatically create this release and utilize the commit history from the `master` branch to generate release notes.
    * `pypi-release` Finally, we utilize the release artifacts created from the GitHub release job of the workflow to publish a new version of Ersilia on the Python Packaging Index (PyPI).
    * The completion of this release workflow triggers the Docker build workflow mentioned above.
*   `Conda-Forge Release`

    Ersilia is released on conda-forge using by way of updating its [feedstock](https://github.com/conda-forge/ersilia-feedstock) maintained within the conda-forge organization on GitHub. Upon the creation of a new GitHub release, an automated pull request is created in this feedstock repository which is then manually reviewed and merged by an Ersilia maintainer.

### Community

* **Creating a Model Request**

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
