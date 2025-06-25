---
description: >-
  Ersilia relies on GitHub Actions for CI/CD. Here is a high-level summary of
  workflows involved in maintaining the Ersilia CLI and code-bases from all the
  models available via the Hub.
---

# CI/CD workflows

## Ersilia CLI

These workflows ensure that we push quality code to [Ersilia](https://github.com/ersilia-os/ersilia), incorporating testing by individual components as well as end-to-end.&#x20;

### Testing

#### Testing Ersilia on a Pull Request/Push

* Objective: Ensure Ersilia continues to work with any new code merging from contributors and maintainers to the master branch. Comprehensive tests include several Python versions, supported platforms and testing several models. This workflow also runs the unit and integration test suite for the repository.
* Workflow name: Test Ersilia codebase
* Workflow file: [tests\_and\_cleanup.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/tests_and_cleanup.yml)
* Runs on: automatically when a PR is open on the Ersilia master branch or a push is made directly to master.
* Jobs:
  * `install-ersilia`, installs Ersilia from Python 3.8 to Python 3.12 on an Ubuntu base.
  * `test-docker`  is the next step to follow after we confirm that Ersilia can be installed in environments created with different Python versions. This step tests whether dockerized models can be fetched using the ersilia CLI, ie the happy flow works without issues.
  * `test-os`  Since Ersilia is supported only on Linux distributions, and MacOS, this job ensures that Ersilia can be installed and run on both of these platforms, presently only utilizing a Python 3.10 environment.&#x20;
  * `run-pytest`  This is the main step that carries out unit testing within Ersilia ensuring that different components continue to work without issues.
  * `run-cli-test-single` This an integration level test that ensures different components within Ersilia interact cohesively in running a single model using the Ersilia CLI. More details on how this is implemented can be found [here](testing-playground.md).
  * `run-cli-test-multiple` Similar to the job described above, this job tests the functionality in Ersilia to simultaneously run multiple models. More details around implementation of the specific tests and the testing frameworks used can be found [here](testing-playground.md).
  * `test-colab-notebook` This test ensures that Ersilia can run in IPython notebooks, such as Jupyter or Google Colab. Specifically, this test checks whether a notebook containing Ersilia code can be executed successfully.
  * `update-model-request-template` This job within the workflow only runs when the PR or the commit on `master` branch change the file containing [model tags](https://github.com/ersilia-os/ersilia/blob/master/ersilia/hub/content/metadata/tag.txt). We expect model tags to be completely sorted because they are used to populate the [Model Request Issue](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&projects=\&template=model_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) (more on this below). The job sorts the tags file if it is changed, and then updates the Model Request Issue [template](https://github.com/ersilia-os/ersilia/blob/master/.github/ISSUE_TEMPLATE/model_request.yml) with the changes, and commits them back to the master branch.

### Packaging

#### Creating Docker image of Ersilia's CLI

* Objective: Package Ersilia into an image to use as a base image to containerise individual models using BentoML.
* Workflow name: Upload Ersilia base image to DockerHub
* Workflow file: [ersilia-base-image-to-dockerhub.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/ersilia-base-image-to-dockerhub.yml)
* Runs on: when code is pushed on to the `master` branch ro when a new release is created from Ersilia at the beginning of each month (more on this below)
* Jobs:
  * `upload_ersilia_base_to_dockerhub` builds two versions of the Ersilia image,  a lighter version tagged as _`ersiliaos/base-multistage-condapack` ,_ and a legacy bloated version tagged as _`ersiliaos/base-legacy-bentoml` ._&#x20;

#### Ersilia Release

* Objective: Distribute Ersilia as source, as a Docker image, as a PyPI package, and as a conda-forge package.&#x20;
* Workflow name: Ersilia Release
* Workflow file: [publish.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/publish.yaml)
* Runs on: At the beginning of each month 3:00 AM UTC
* Jobs:
  * `version` increments the version of the codebase and creates a new release tag, while also updating this information in the [CITATION.cff](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff), and [CODEMETA.json](https://github.com/ersilia-os/ersilia/blob/master/codemeta.json) files.&#x20;
  * `gh-release` Utilizing the tag created in the previous job, a GitHub release is created. The GitHub CLI is used to automatically create this release and utilize the commit history from the `master` branch to generate release notes.
  * `pypi-release` Finally, we utilize the release artifacts created from the GitHub release job of the workflow to publish a new version of Ersilia on the Python Packaging Index (PyPI).
  * The completion of this release workflow triggers the Docker build workflow mentioned above.
  * `Conda-Forge Release` Ersilia is released on conda-forge using by way of updating its [feedstock](https://github.com/conda-forge/ersilia-feedstock) maintained within the conda-forge organization on GitHub. Upon the creation of a new GitHub release, an automated pull request is created in this feedstock repository which is then manually reviewed and merged by an Ersilia maintainer.

### Community

#### New Model Request

* Objective: Contribute new models to the Ersilia Model Hub
* Workflow name: Approve Command Dispatch
* Workflow file: [approve-dispatch.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/approve-dispatch.yml)
* Runs on: When an Ersilia maintainer comments `/approve`  on a [Model Request](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&projects=\&template=model_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) issue submitted by a model contributor. A small action named `approve-trigger.yml` first checks if we are effectively on a Model request issue.
* Jobs:&#x20;
  * `approve-command-dispatch`  The Model Request issue is parsed to ensure all required fields are filled in. If the metadata is correct, a new model repository is created from the [eos-template](https://github.com/ersilia-os/eos-template) and a record is generated in Airtable. The repository includes all the files in eos-template and a mock Git LFS object to allow for uploading Git LFS objects in the repository if required. The issue creator is notified to prompt them to start contributing to the model. Read on more about model contributions in our Model Contribution [guidelines](../model-contribution/example-of-the-model-incorporation-workflow.md).&#x20;

## Ersilia Models

Ersilia models include a series of workflow callers for the reusable workflows available in the [ersilia-model-workflows](https://github.com/ersilia-os/ersilia-model-workflows) repository:

### Test and upload model

* Objective: Test a model using its source code and, if correct, update it across all the platforms (GitHub, S3 and Dockerhub)
* Workflow file: [upload-model.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/upload-model.yml)
* Runs on: whenever code is pushed to the main branch of the repository, either directly by a maintainer or by way of merging a pull request. This workflow can also be triggered manually by a maintainer.
* Jobs:
  * `test-model-source`: uses the reusable workflow of the same name, which installs Ersilia and tests the source code of the model using the `test --shallow` command. It captures and updates the relevant metadata to Airtable
  * `upload-model-to-s3`: uses the reusable workflow of the same name and uploads the model to S3 if the test of the source code has passes
  * `upload-ersilia-pack`: uses the reusable workflow of the same name and is triggered if the upload to S3 is successful. A Docker image is built for both AMD64 and ARM64 architectures and tagged as dev. By default we first try to build a FastAPI packaged model
  * `upload-bentoml-multistage`: uses the upload-bentoml.yml workflow with the version multistage-condapack if the build with ersilia-pack has failed. It also tries to build a Docker image for both AMD64 and ARM64 architectures using the ersilia-maintained BentoML.
  * `upload-bentoml-legacy`: uses the upload-bentoml.yml workflow with the version legacy if the build with multistage-condapack has failed. It also tries to build a Docker image for both AMD64 and ARM64 architectures using the legacy version of BentoML.

### Test model image

* Objective: test the dev version of a model image
* Workflow file: [test-model-image.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/test-model-image.yml)
* Runs on: whenever the test and upload model has successfully completed one of the three Docker build options
* Jobs:
  * `test-image`: using the reusable workflow test-model-image, it tests the dev versions of the image for both architectures if available using the `test --deep` command. If successful, retags the image(s) as latest.
  * `post-upload`: the post-model-upload workflow is triggered if the image testing is successful. It reads the output of the test command (.json) and captures the relevant metadata (image size, computational performance etc). It also changes the status of the model to Ready. Finally it updates Airtable and the metadata and readme files in the model repository.

### Test model PR

* Objective: test that the model incorporation is complete
* Workflow file: [test-model-pr.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/test-model-pr.yml)
* Runs on: whenever a PR is opened to the model repository
* Jobs:&#x20;
  * `test-model-pr`: using the workflow of the same name, the workflow installs ersilia, clones the code in the PR and performs a `test --shallow` to check the model is able to work within Ersilia.&#x20;

## Production and serving

* Ersilia Statistics&#x20;
* Ersilia Maintenance
