# Test command

The Ersilia test command is designed to automatize model testing and validation at incorporation time and during routine maintenance of the models.&#x20;

## TL:DR

The test command is a CLI command on the Ersilia Model Hub that automatically performs several checks on an individual model.

To run the test command you need to install Ersilia with the extra packages required for testing:

```bash
conda activate ersilia
pip install ersilia[test]
ersilia test model_id
```

The test command has three levels of complexity:

* Surface: high level tests to ensure all the necessary files are available and metadata is compliant with Ersilia standards. It does not actually run the model, and is designed to be a quick maintenance for models.
* Shallow: will run all the tests from the basic command and then test the model through Ersilia's CLI. Designed to be run by model contributors prior to incorporating the model and when any change is introduced in the model.
* Deep: performs all shallow tests and in addition calculates performance metrics for the model. Designed to be used only by Ersilia maintainers when a new model is incorporated or significant changes are introduced.

Hence, a test command could look like:

```
ersilia test model_id --surface
ersilia test model_id --shallow
ersilia test model_id --deep
```

{% hint style="info" %}
In addition to the printed output in the terminal, the test command produces a .json report stored in the directory where the command is run under the name `<model_id>-test.json`
{% endhint %}

## Surface test

### Usage

The model is downloaded (not fetched) from its online storage in S3 or Github. By default, if no flag is specified the model will be fetched from GitHub. In addition, the basic tests can also be performed from a local directory, for example when a contributor is incorporating a new model.

```
ersilia test model_id [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket)</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Uses a model stored locally on the indicated path.</td></tr></tbody></table>

### Tests performed

* Model information: the model metadata is available in the `metadata.json` or `metadata.yml` and in the correct format. If the model is not yet fully incorporated in Ersilia, fields like S3 URL or Docker Architecture will not exist. Importantly those are marked as Not Passed, instead of Failed.
* Model files: all required model files exist for either of Ersilia package modes:
  * BentoML packaging: Dockerfile, metadata.json, run.sh, service.py, pack.py, README.md, LICENSE.
  * FastAPI packaging: install.yml, metadata.yml. input.csv, output.csv, run.sh, README.md, LICENSE
* Example model files:
  * run\_input.csv and run\_output.csv
  * run\_columns.csv: output columns and its information are properly specified
* Model directory size: calculates the size of the model directory (including model checkpoints)
* Dependencies: ensures all dependencies specified are pinned to a specific version.

### Outputs

The terminal will print four tables, one per each type of test specified above, and whether each check has PASSED or FAILED. In the `.json` file, the tests appear as True (passed) or False (failed).

## Shallow test

### Usage

This command will not only download the model from the specified source but also test that it actually works by fetching, serving and running it via Ersilia's CLI.&#x20;

```
ersilia test model_id --shallow [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation and then fetches it from the created folder.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket) and then fetches it from the created folder.</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Fetches a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Fetches the model from DockerHub, and in parallel downloads it from GitHub to perform the basic tests</td></tr></tbody></table>

{% hint style="info" %}
To perform the basic test checks the model needs to be Downloaded, not Fetched, so we can calculate properties like Model Directory Size. Hence, the --shallow and --deep tests use a combination of downloading and fetching.
{% endhint %}

### Tests performed

In addition to the basic tests:

* Environment size: total size of the enviroment created when installing the model dependencies. Only available if fetching --from\_dir/from\_s3/from\_github.
* Container size: total size of the Docker Image downloaded from DockerHub. Only available if fetching --from\_dockerhub.
* Output correctness: the test module performs several runs and ensures:
  * The output contains the columns defined in run\_columns.csv
  * The output can be passed in all accepted formats (single molecule, list or `.csv`) saved in all available formats (`.csv`, `.json`, `.h5`)
  * The output is consistent between runs (for the same molecule, same result or small divergence in non-stochastic models)
  * The output is consistent between running the model via Ersilia or directly from the run-sh command.
  * The output does not contain a majority of null values.

### Outputs

The basic checks are printed in the terminal first and then once the basic tests are completed the model is fetched, served and run, providing the following tables:

* Validation and size check results: will evaluate the model size (environment or image) and the consistency of the model output. If an example input file is not included in the model, it will be created from a random sampling of Ersilia Maintained Inputs. The consistency will be calculated for both string and numerical outputs using scores like `rmse` and `spearmanr` (Spearman correlation coefficient with associated p-value).
* Consistency summary between Ersilia and Bash execution options: will compare the consistency of running the model via Ersilia or directly running the run.sh file on a virtual environment.
* Model output content validation summary: ensures all the combinations between input and output are valid and working.

## Deep test

### Model source

The command will perform the basic and shallow tests and in addition run predictions for one, fifty and one hundred inputs and calculate the computational performance.

```
ersilia test model_id --deep [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation and then fetches it from the created folder.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket) and then fetches it from the created folder.</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Fetches a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Fetches the model from DockerHub, and in parallel downloads it from GitHub to perform the basic tests</td></tr></tbody></table>

### Tests performed

Computational performance assessment: after serving the models, they get executed for 1, 50, and 100 inputs. Model performance (seconds) will be recorded and reported using `wall clock` .

### Outputs

The same tables as in the basic and shallow test and an additional model performance table evaluating the performance (time (s)) taken to run inputs of lenght:

* 1
* 10
* 100
* 1000
* 10000

## Detailed Methods

The mechanism involves several services and classes working together to ensure the model's functionality and reliability. The main components are:

* `RunnerService`: Manages the execution of model tests and checks.
* `InspectService`: Inspects models and their configurations.
* `CheckService`: Performs various high level checks mentioned above checks on the model.
* `IOService`: Handles input/output operations related to model testing.
* `ModelTester`: A high-level class that orchestrates the testing process.

The process typically involves:

1. Setting up the environment and fetching the model repository.
2. Running various checks to ensure the model's integrity.
3. Generating tables and logs of the results.
4. Cleaning up temporary files and directories (model directotry only if the `--remove` flag is enabled).

***
