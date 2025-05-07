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

* Basic: high level tests to ensure all the necessary files are available and metadata is compliant with Ersilia standards. It does not actually run the model, and is designed to be a quick maintenance for models
* Surface: performs all basic tests and a simple run of the model fetched from the specified source through Ersilia's CLI. Designed to be run as a quick check during model maintenance.
* Shallow: will run all the tests from the surface command and then test the model more thoroughly by testing all input types (`string`, `list`, `csv`) and output types (`csv`, `json`, `h5`). Also ensures the consistency of results between runs. Designed to be run by model contributors prior to incorporating the model and when any change is introduced in the model.
* Deep: performs all shallow tests and in addition calculates performance metrics for the model. Designed to be used only by Ersilia maintainers when a new model is incorporated or significant changes are introduced.

Hence, a test command could look like:

```bash
ersilia test model_id --inspect #Basic run
ersilia test model_id --surface #Surface check
ersilia test model_id --shallow #Shallow check
ersilia test model_id --deep #Deep check
```

{% hint style="info" %}
In addition to the printed output in the terminal, the test command produces a .json report stored in the directory where the command is run under the name `<model_id>-test.json`
{% endhint %}

## Inspect test

### Usage

The model is downloaded (not fetched) from its online storage in S3 or Github. By default, if no flag is specified the model will be fetched from GitHub. In addition, the basic tests can also be performed from a local directory, for example when a contributor is incorporating a new model.

```
ersilia test model_id --inspect [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket)</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Uses a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>It will default back to from_github, as the model needs to be downloaded</td></tr></tbody></table>

### Tests performed

* Metadata checks: the model metadata is available in the `metadata.json` or `metadata.yml` and in the correct format. If the model is not yet fully incorporated in Ersilia, fields like S3 URL or Docker Architecture will not exist. Importantly those are marked as Not Passed, instead of Failed.
* Model file checks: all required model files exist for either of Ersilia package modes:
  * BentoML packaging: Dockerfile, metadata.json, run.sh, service.py, pack.py, README.md, LICENSE.
  * FastAPI packaging: install.yml, metadata.yml, README.md, LICENSE files, model/framework/examples/run\_input.csv, model/framework/examples/run\_output.csv, model/framework/columns/run\_columns.csv, model/framework/run.sh,&#x20;
* File validity check: checks that the following files have the expected structure:
  * Columns: the columns specified in the `run_columns.csv` and the `run_output.csv` coincide.
  * Dockerfile/Install\_yml: the dockerfile or install.yml contains the proper formatting of dependencies (see more in the [Model Contribution](../model-contribution/) section)
* Model directory size: calculates the size of the model directory (including model checkpoints)

### Outputs

The terminal will print four tables, one per each type of test specified above, and whether each check has PASSED or FAILED. In the `.json` file, the tests appear as True (passed) or False (failed). The `-v` flag can always be used to see more information on the terminal.

## Surface test

### Usage

The model is downloaded (not fetched) from its online storage in S3 or Github, and fetched via the specified source (including Dockerhub). By default, if no flag is specified the model will be fetched from GitHub. In addition, the basic tests can also be performed from a local directory, for example when a contributor is incorporating a new model.

```
ersilia test model_id --surface [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket)</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Uses a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Downloads the model from GitHub and fetches it via DockerHub for the simple run</td></tr></tbody></table>

### Tests performed

All the basic tests and in addition

* Model size check: it also calculates the environment size (if fetched `--from_github`, `--from_s3` or `--from_dir`) or image size (if fetched `--from_dockerhub`)
* Model run check: fetches the model through Ersilia's CLI from the specified source and then runs the `run_input.csv`. Outputs the result in `.csv` format and ensures that the result is not all None's and that the columns match the results in `run_output.csv`

### Outputs

The terminal will print five tables, one per each type of test specified above, and whether each check has PASSED or FAILED. In the `.json` file, the tests appear as True (passed) or False (failed). The `-v` flag can always be used to see more information on the terminal.

## Shallow test

### Usage

Similar to the surface but performs more tests once the model is fetched.

```
ersilia test model_id --shallow [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation and then fetches it from the created folder.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket) and then fetches it from the created folder.</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Fetches a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Fetches the model from DockerHub, and in parallel downloads it from GitHub to perform the basic tests</td></tr></tbody></table>

{% hint style="info" %}
To perform the basic test checks the model needs to be Downloaded, not Fetched, so we can calculate properties like Model Directory Size. Hence, the --shallow and --deep tests use a combination of downloading and fetching.
{% endhint %}

### Tests performed

In addition to the basic tests:

* Input output check: the test module performs several runs and ensures that the output can be passed in all accepted formats (single molecule, list or `.csv`) saved in all available formats (`.csv`, `.json`, `.h5`)
* Model output consistency check: the output is consistent between runs (for the same molecule, same result or small divergence in non-stochastic models). The consistency will be calculated for both string and numerical outputs using scores like `rmse` and `spearmanr` (Spearman correlation coefficient with associated p-value). The `example/run_input.csv` file is used for this check.
* Consistency summary between ersilia and bash execution: the output is consistent between running the model via Ersilia or directly from the `run.sh` bash file in `model/framework`. The `example/run_input.csv` file is used for this check.

### Outputs

The terminal will print eight tables, one per each type of test specified above, and whether each check has PASSED or FAILED. In the `.json` file, the tests appear as True (passed) or False (failed). The `-v` flag can always be used to see more information on the terminal.

## Deep test

### Model source

The command will perform the basic and shallow tests and in addition run predictions for one, fifty and one hundred inputs and calculate the computational performance.

```
ersilia test model_id --deep [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation and then fetches it from the created folder.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket) and then fetches it from the created folder.</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Fetches a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Fetches the model from DockerHub, and in parallel downloads it from GitHub to perform the basic tests</td></tr></tbody></table>

### Tests performed

* Computational performance assessment: after serving the models, they get executed for 1, 10, 100, 1000 and 10000 inputs. Model performance (seconds) will be recorded and reported using `wall clock` . By default, a `deterministic` flag is enabled for the `example` command, ensuring all models use the same list of molecules for the test command.

### Outputs

The same tables as in the basic, surface and shallow test and an additional model performance table evaluating the performance (time (s)) taken to run inputs of each length.

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
4. Cleaning up temporary files and directories (model directory only if the `--remove` flag is enabled).

***
