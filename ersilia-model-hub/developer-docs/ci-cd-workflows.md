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

#### Testing Ersilia on a Pull Request

* Objective: Ensure Ersilia continues to work with any new code merging from contributors and maintainers to the master branch
* Workflow name: Deploy and test ersilia on PR
* Workflow file: [pr\_check.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/pr_check.yml)
* Runs on: automatically when a PR is open on the Ersilia master branch.
* Jobs:&#x20;
  * `Build`: creates a conda environment (py3.10 currently) and installs the Ersilia codebase from the open PR. Then it tests the following functions: catalog and model running (fetch, serve, info and run) from GitHub, S3 and DockerHub. A single molecule is passed as input, and the output is analysed.

#### Testing and cleaning up

* Objective: Ensure that Ersilia is installable across supported Python environments, supported platforms, both through source and while interacting with Dockerized models. This workflow also runs the unit and integration test suite for the repository.
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
* Run: when a PR is opened against the master branch, and again when code is pushed to the master branch.

### Packaging

#### Creating Docker image of Ersilia's CLI

* Objective: Package Ersilia into an image to use as a base image to containerise individual models as well as in development containers.
* Workflow name: Upload Ersilia base image to DockerHub
* Workflow file: [ersilia-base-image-to-dockerhub.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/ersilia-base-image-to-dockerhub.yml)
* Jobs:
  * `upload_ersilia_base_to_dockerhub` builds two versions of the Ersilia image,  a lighter version tagged as _`ersiliaos/base-multistage-condapack` ,_ and a legacy much bloated version tagged as _`ersiliaos/base-legacy-bentoml` ._&#x20;
* Run: configured to trigger in multiple ways - whenever code is pushed on to the `master` branch, when a new release is created from Ersilia at the beginning of each month (more on this below), and also manually by developers and maintainers.&#x20;

#### Ersilia Release

* Objective: Distribute Ersilia as source, as a Docker image, as a PyPI package, and as a conda-forge package.&#x20;
* Workflow name: Ersilia Release
* Workflow file: [publish.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/publish.yaml)
* Jobs:
  * `version` This is the job that increments the version of the code-base, and creates a new release tag, while also updating this information in the [CITATION.cff](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff), and [CODEMETA.json](https://github.com/ersilia-os/ersilia/blob/master/codemeta.json) files.&#x20;
  * `gh-release` Utilizing the tag created in the previous job, a GitHub release is created in this job of the workflow. The GitHub CLI is used to automatically create this release and utilize the commit history from the `master` branch to generate release notes.
  * `pypi-release` Finally, we utilize the release artifacts created from the GitHub release job of the workflow to publish a new version of Ersilia on the Python Packaging Index (PyPI).
  * The completion of this release workflow triggers the Docker build workflow mentioned above.
  *   `Conda-Forge Release`

      Ersilia is released on conda-forge using by way of updating its [feedstock](https://github.com/conda-forge/ersilia-feedstock) maintained within the conda-forge organization on GitHub. Upon the creation of a new GitHub release, an automated pull request is created in this feedstock repository which is then manually reviewed and merged by an Ersilia maintainer.
* Run on: At the beginning of each month, a release [workflow](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/publish.yaml) is executed at 3:00 AM UTC publishing a new release for Ersilia using the code present on the `master` branch.&#x20;

### Community

#### New Model Request

* Objective: Contribute new models to the Ersilia Model Hub
* Workflow name: Approve Command Dispatch
* Workflow file: [approve-dispatch.yml](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/approve-dispatch.yml)
* Jobs:&#x20;
  * `approve-command-dispatch`  A model contributor is expected to fill out several fields in the [Model Request](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&projects=\&template=model_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) issue form such as the model's name, its slug, relevant tags, its publication and source code, etc. Once approved, the workflow creates an model repository from the [model-template](https://github.com/ersilia-os/eos-template) utilizing information parsed from the issue body to name the repository, track information about the requested model within Ersilia's metadata stores such as Airtable, and read-only JSON files in S3. A mock LFS object is also instantiated and uploaded to the repository to facilitate model contributors and utilizing Ersilia's [Git LFS quota](https://github.com/homuler/MediaPipeUnityPlugin/issues/475) while adding checking model weights in the repository. Once these steps complete successfully, the original issue creator is notified that a repository has been successfully created for them to start contributing the model. Read on more about model contributions in our Model Contribution [guidelines](../model-contribution/example-of-the-model-incorporation-workflow.md).&#x20;
* Run: When an Ersilia maintainer comments `/approve`  on a [Model Request](https://github.com/ersilia-os/ersilia/issues/new?assignees=\&labels=new-model\&projects=\&template=model_request.yml\&title=%F0%9F%A6%A0+Model+Request%3A+%3Cname%3E) issue submitted by a model contributor.&#x20;

## Model incorporation

### Testing Model Source

#### Directly Modifying the Model Repository

* Objective: Test a model using its source code.
* Workflow name: Model test on push
* Workflow file: [test-model.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/test-model.yml)
* Jobs:
  * `test`  This workflow is responsible for testing the model source and verifying that the model produces outputs and those outputs are functionally correct. The model is first fetched using its local path on the runner, and is then run against a set of example inputs. The model output is checked for [null values](https://github.com/ersilia-os/eos-template/blob/main/.github/scripts/verify_model_outcome.py) as a proxy for making sure that the model actually works.  Thereafter,  the `ersilia test`  command is utilized to perform sanity checks like completeness of model dependencies, metadata, and necessary files, including testing for functional correctness of the model's output. The job is also responsible for updating the model's metadata in Ersilia's metadata stores such as Airtable and read-only JSON files in S3.
* Run: Whenever code is pushed to the main branch of the repository, either directly by a maintainer or by way of merging a pull request. This workflow can also be triggered manually by a maintainer.

#### Creating a Pull Request

* Objective: Test a model using its source code.
* Workflow name: Model Test on PR
* Workflow file: [test-model-pr.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/test-model-pr.yml)
* Jobs:
  * `test` This workflow is quite similar to the one above. It tests a model by fetching it locally from source, and verifying that the model both produces outputs and that they are functionally correct. However, this workflow does not sync any metadata with Ersilia's metadata stores.
* Run: Upon Pull Request creation by the model contributor from the forked model repository.

### Packaging Model Images

{% hint style="info" %}
Docker images for _most_ models within the hub are generally supported on both AMD and ARM chips.
{% endhint %}

#### Upload Model to S3

* Objective: Upload the model source code as a zip file to Ersilia's S3 storage in AWS cloud.
* Workflow name: Upload model to S3
* Workflow file: [upload-model-to-s3.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/upload-model-to-s3.yml)
* Jobs:
  * `upload_model_to_s3`  This workflow is responsible for zipping the model source and [uploading](https://github.com/ersilia-os/ersilia/blob/master/.github/scripts/upload_model_to_s3.py) it to a S3 bucket in Ersilia's AWS cloud. The referenced script performs many actions such as removing git artifacts from the model source directory before preparing it for upload. The workflow also carries out updating metadata for the model in Ersilia's metadata stores, particularly with the URL of the AWS S3 object that gets created by running this workflow.
* Run: Only on the main repository, upon successful completion of the `Model test on push`  workflow.

#### Upload Model to DockerHub

* Objective: Upload a functioning dockerized image for the model
* Workflow name: Upload model to DockerHub
* Workflow file: [upload-model-to-dockerhub.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/upload-model-to-dockerhub.yml)
* Jobs:
  * `upload-ersilia-pack`  This job uses the callable [`upload-ersilia-pack.yml`](ci-cd-workflows.md#upload-ersilia-pack) workflow to build and upload the Ersilia Pack Dockerized model. The called workflow inherits any necessary secrets from the repository.
  * `upload-multi-stage-condapack`  This job runs if the `upload-ersilia-pack` job fails. It uses the upload-bentoml.yml workflow to build and upload a multi-stage CondaPack Dockerized model. Because it depends on the success or failure status of the `upload-ersilia-pack`  it cannot start before the previous job is finished. This job uses the callable workflow [`upload-bentoml.yml`](ci-cd-workflows.md#upload-bentoml) with the input "multistage-condapack". More on how models are Dockerized within Ersilia can be found [here](dockerization-of-ersilia-models.md).&#x20;
  * `upload-legacy-bentoml`  This job runs if the `upload-multi-stage-condapack` job fails. It uses the callable [`upload-bentoml.yml`](ci-cd-workflows.md#upload-bentoml) workflow to build and upload a legacy BentoML Dockerized model. More information on this terminology (and some history) can be found [here](dockerization-of-ersilia-models.md).
* Run: When the `Upload Model to S3`  workflow is successfully completed.

#### Upload Ersilia Pack

* Objective: This workflow is called by [another workflow](ci-cd-workflows.md#upload-model-to-dockerhub) to build and upload the Ersilia Pack Dockerized model
* Workflow name: Upload Ersilia Pack Dockerized Model
* Workflow file: [upload-ersilia-pack.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/upload-ersilia-pack.yml)
* Jobs:&#x20;
  * `build-ersilia-pack-image` The job fetches the [correct Dockerfile](https://github.com/ersilia-os/eos-template/blob/main/.github/scripts/resolve_dockerfile.py) based on the model's Python version requirement, and uses that Dockerfile to build an image for the model. Thereafter, the model is fetched using its locally built image, served with tracking enabled, and run against a set of example inputs. The output is checked to see if it contains any [not null](https://github.com/ersilia-os/eos-template/blob/main/.github/scripts/verify_model_outcome.py) values as a proxy measure for testing whether the model image got built successfully or not. Note that we only build an AMD image for testing, and if the test passes, we build the ARM image as well, and push them both to DockerHub. Note that if the ARM image fails to build, we store that in `arch.txt`  build artifact.&#x20;
* Run: When called by `Upload model to DockerHub`  workflow.

#### Upload BentoML&#x20;

* Objective: This workflow is called by other workflows and requires an input version to determine the type of BentoML Dockerized model to build.
* Workflow name: Upload BentoML Dockerized Model with or without Multi-stage Conda Pack
* Workflow file: [upload-bentoml.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/upload-bentoml.yml)
* Jobs:
  * `build-bentoml-image` This job is only run if the model repository structure supports building a model using BentoML, ie, it should not contain an `install.yml`  file. More on why this happens can be found [here](dockerization-of-ersilia-models.md). This workflow is quite similar to the workflow above, in that it [generates](https://github.com/ersilia-os/ersilia/blob/master/.github/scripts/place_a_dockerfile_in_current_eos_repo.py) the Dockerfile for the model repository, and builds and tests the model image. First only the AMD image is built, which is then tested for [not null](https://github.com/ersilia-os/eos-template/blob/main/.github/scripts/verify_model_outcome.py) outcomes as a proxy measure for a functional model image, and if this test passes, the ARM image also gets built. After successful building of both versions, the images get pushed to DockerHub. Note that if the ARM image fails to build, we store that in `arch.txt`  build artifact.&#x20;
* Run: When called by `Upload model to DockerHub`  workflow.

#### Post Model Upload Actions

* Objectives: Update metadata for the model, and notify users to test the model.
* Workflow name: Post Model Upload actions
* Workflow file: [post-model-upload.yml](https://github.com/ersilia-os/eos-template/blob/main/.github/workflows/post-model-upload.yml)
* Jobs:
  * `post-model-upload`  This job utilizes the `arch.txt` build artifact from the `Upload Model to DockerHub`  job to update the model metadata file in the repository signifying the image architectures supported by the model, as well as updates that information in Ersilia's metadata stores such as Airtable, and read-only JSON files in S3. Thereafter, a test issue is created in the model's repository notifying an Ersilia maintainer to test the model. If the Ersilia maintainer finds the model to be working, then at this stage we conclude the model to have been successfully incorporated.
* Run: Upon completion of the `Upload model to DockerHub`  workflow.

## Production and serving

* Ersilia Statistics&#x20;
* Ersilia Maintenance
