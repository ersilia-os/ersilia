# Model Tester

This functionality provides tools and tests for running and validating individual AI/ML models that are under development locally or even deployed remotely. Models in Ersilia are packed with two methods, . The `test` command supports AI/ML models packaged with both `bentoml` or `ersilia-pack` methods (see <mark style="color:red;">this page</mark> for more information).

## TL:DR

The test command is a CLI command on the Ersilia Model Hub that automatically performs several checks on an individual model.

To run the test command you need to install Ersilia with the extra packages required for testing:

```bash
conda activate ersilia
pip install ersilia[test]
ersilia test model_id
```

The test command has three levels of complexity:

* Basic: high level tests to ensure all the necessary files are available and metadata is compliant with Ersilia standards. It does not actually run the model, and is designed to be a quick maintenance for models.
* Shallow: will run all the tests from the basic command and then test the model through Ersilia's CLI. Designed to be run by model contributors prior to incorporating the model and when any change is introduced in the model.
* Deep: performs all shallow tests and in addition calculates performance metrics for the model. Designed to be used only by Ersilia maintainers when a new model is incorporated or significant changes are introduced.

Hence, a test command could look like:

```
ersilia test model_id
ersilia test model_id --shallow
ersilia test model_id --deep
```

## Basic test

### Usage

The model is downloaded (not fetched) from its online storage in S3 or Github. By default, if no flag is specified the model will be fetched from GitHub. In addition, the basic tests can also be performed from a local directory, for example when a contributor is incorporating a new model.

```
ersilia test model_id [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket)</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Uses a model stored locally on the indicated path.</td></tr><tr><td>--as_json</td><td>False</td><td>Saves the output of the tests performed in a machine-readable json file.</td></tr></tbody></table>

### Tests performed

* Model information: the model metadata is available in the `metadata.json` or `metadata.yml` and in the correct format. If the model is not yet fully incorporated in Ersilia, fields like S3 URL or Docker Architecture will not exist.
* Model files: all required model files exist for either of Ersilia package modes:
  * BentoML packaging: Dockerfile, metadata.json, run.sh, service.py, pack.py, README.md, LICENSE.
  * FastAPI packaging: install.yml, metadata.yml. input.csv, output.csv, run.sh, README.md, LICENSE
* Model directory size: calculates the size of the model directory (including model checkpoints)
* Dependencies: ensures all dependencies specified are pinned to a specific version.

### Outputs

The terminal will print four tables, one per each type of test specified above, and whether each check has PASSED or FAILED.&#x20;

If the `--as_json` flag is passed, the output will be printed in the terminal and also saved in the directory where you are at with the name `model_id-test.json`. In the Json file, the tests appear as True (passed) or False (failed).

## Shallow test

### Usage

This command will not only download the model from the specified source but also test that it actually works by fetching, serving and running it via Ersilia's CLI.&#x20;

```
ersilia test model_id --shallow [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation and then fetches it from the created folder.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket) and then fetches it from the created folder.</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Fetches a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Fetches the model from DockerHub, and in parallel downloads it from GitHub to perform the basic tests</td></tr><tr><td>as_json</td><td>False</td><td>Saves the output of the tests performed in a machine-readable json file.</td></tr></tbody></table>

{% hint style="info" %}
To perform the basic test checks the model needs to be Downloaded, not Fetched, so we can calculate properties like Model Directory Size. Hence, the --shallow and --deep tests use a combination of downloading and fetching.
{% endhint %}

### Tests performed

In addition to the basic tests:

* Environment size: total size of the enviroment created when installing the model dependencies. Only available if fetching --from\_dir/from\_s3/from\_github.
* Container size: total size of the Docker Image downloaded from DockerHub. Only available if fetching --from\_dockerhub.
* Output correctness: the test module performs several runs and ensures:
  * The output can be passed in all accepted formats (single molecule, list or `.csv`) saved in all available formats (`.csv`, `.json`, `.h5`)
  * The output is consistent between runs (for the same molecule, same result or small divergence in non-stochastic models)
  * The output does not contain a majority of null values.

### Outputs

The basic checks are printed in the terminal first and then once the basic tests are completed the model is fetched, served and run, providing the following tables:

* Validation and size check results: will evaluate the model size (environment or image) and the consistency of the model output. If an example input file is not included in the model, it will be created from a random sampling of Ersilia Maintained Inputs. The consistency will be calculated for both string and numerical outputs using scores like `rmse` and `spearmanr` (Spearman correlation coefficient with associated p-value).
* Consistency summary between Ersilia and Bash execution options: will compare the consistency of running the model via Ersilia or directly running the run.sh file on a virtual environment.
* Model output content validation summary: ensures all the combinations between input and output are valid and working.

If the `--as_json` flag is passed, the output will be printed in the terminal and also saved in the directory where you are at with the name `model_id-test.json`. In the Json file, the tests appear as True (passed) or False (failed).

## Deep test

### Model source

The command will perform the basic and shallow tests and in addition run predictions for one, fifty and one hundred inputs and calculate the computational performance.

```
ersilia test model_id --deep [OPT?]
```

<table><thead><tr><th width="220">Flag</th><th width="129">Default</th><th>Description</th></tr></thead><tbody><tr><td>--from_github</td><td>True</td><td>Downloads the model from its repository on the ersilia-os organisation and then fetches it from the created folder.</td></tr><tr><td>--from_s3</td><td>False</td><td>Downloads the model from its storage on the cloud (S3 Bucket) and then fetches it from the created folder.</td></tr><tr><td>--from_dir [path/to/dir]</td><td>False</td><td>Fetches a model stored locally on the indicated path.</td></tr><tr><td>--from_dockerhub</td><td>False</td><td>Fetches the model from DockerHub, and in parallel downloads it from GitHub to perform the basic tests</td></tr><tr><td>as_json</td><td>False</td><td>Saves the output of the tests performed in a machine-readable json file.</td></tr></tbody></table>

### Tests performed

Computational performance assessment: after serving the models, they get executed for 1, 50, and 100 inputs. Model performance (seconds) will be recorded and reported using `wall clock` .

### Outputs

The same tables as in the basic and shallow test and an additional model performance table.

If the `--as_json` flag is passed, the output will be printed in the terminal and also saved in the directory where you are at with the name `model_id-test.json`. In the Json file, the tests appear as True (passed) or False (failed).

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

### Key Methods

#### `check_single_input(output_file, run_model, run_example)`

* **Purpose**: Validates a single model input against provided configurations.
* **How it works**:
  * Generates a single input sample using `run_example`.
  * Runs the model with the generated input using `run_model`.
  * Reads the output CSV file to ensure it contains valid results.
*   **Usage**:

    ```python
    check_service.check_single_input(output_file, run_model, run_example)
    ```

#### `_generate_table(title, headers, rows)`

* **Purpose**: Generates a table with the given title, headers, and rows.
* **How it works**:
  * Uses `rich.Table` to create a formatted table.
  * Adds columns and rows based on the provided headers and data.
  * Prints the table to the console.
*   **Usage**:

    ```python
    ios_service._generate_table("Title", ["Header1", "Header2"], [["Row1Col1", "Row1Col2"], ["Row2Col1", "Row2Col2"]])
    ```

#### `_clear_folders()`

* **Purpose**: Removes directories if the `remove` flag is set.
* **How it works**:
  * Uses `SetupService.run_command()` to execute Linux commands for folder cleanup.
  * Deletes the specified directories.
*   **Usage**:

    ```python
    runner_service._clear_folders()
    ```

#### `check_example_input(output_file, run_model, run_example)`

* **Purpose**: Ensures the model processes example data consistently.
* **How it works**:
  * Generates example input samples using `run_example`.
  * Runs the model with the generated input using `run_model`.
  * Verifies that the output is consistent and valid.
*   **Usage**:

    ```python
    check_service.check_example_input(output_file, run_model, run_example)
    ```

#### `check_consistent_output(run_example, run_model)`

* **Purpose**: Verifies that outputs match expected formats or zero differences.
* **How it works**:
  * Generates example input samples using `run_example`.
  * Runs the model twice with the same input using `run_model`.
  * Compares the two outputs to ensure consistency.
*   **Usage**:

    ```python
    check_service.check_consistent_output(run_example, run_model)
    ```

#### `run_bash()`

* **Purpose**: Invokes shell commands to finalize or clean up tests.
* **How it works**:
  * Runs a bash script to execute model tests.
  * Compares the outputs generated by the bash script and the model.
*   **Usage**:

    ```python
    runner_service.run_bash()
    ```

***

### Class: `ModelTester`

Inherits from `ErsiliaBase`.\
**Parameters:**

* `model_id`: Identifies which model is tested.
* `level`: Controls verbosity or testing depth.
* `dir`: Directory for storing temporary files.
* `inspect`, `remote`, `remove`: Booleans that configure optional behavior.

**Usage:**

1. Create an instance with appropriate parameters.
2. Run the desired methods (e.g., checking single input, generating output tables).
3. Review logs or console output for results.

**Example**:

```python
tester = ModelTester(model_id="eosxxxx", level="deep", dir="/path/to/dir", inspect=True, remote=False, remove=True)
tester.setup()
tester.run(output_file="result.csv")
```

***

### Class: `RunnerService`

Manages the execution of model tests and checks.

**Parameters:**

* `model_id`: Identifier of the model.
* `logger`: Logger for logging messages.
* `ios_service`: Instance of `IOService` for handling input/output operations.
* `checkup_service`: Instance of `CheckService` for performing various checks on the model.
* `setup_service`: Instance of `SetupService` for setting up the environment and fetching the model repository.
* `model_path`: Path to the model.
* `level`: Level of checks to perform.
* `dir`: Directory where the model repository is located.
* `remote`: Flag indicating whether to fetch the repository from a remote source.
* `inspect`: Flag indicating whether to perform inspection checks.
* `remove`: Flag indicating whether to remove the model directory after tests.
* `inspecter`: Instance of `InspectService` for inspecting models and their configurations.

**Methods**:

* `run_model(input, output, batch)`: Runs the model with the given input and output parameters.
* `fetch()`: Fetches the model repository from the specified directory.
* `run_exampe(n_samples, file_name, simple, try_predefined)`: Generates example input samples for the model.
* `run_bash()`: Runs the model using a bash script and compares the outputs for consistency.
* `make_output(elapsed_time)`: Generates the final output table with the test results.
* `run(output_file)`: Runs the model tests and checks.

**Example**:

```python
runner = RunnerService(model_id="eosxxxx", logger=logger, ios_service=ios, checkup_service=check_service, setup_service=setup_service, model_path="/path/to/model", level="deep", dir="/path/to/dir", remote=False, inspect=True, remove=True, inspecter=inspect_service)
runner.run(output_file="result.csv")
```

***

### Class: `InspectService`

Service for inspecting models and their configurations.

**Parameters:**

* `dir`: Directory where the model is located.
* `model`: Model identifier.
* `config_json`: Path to the configuration JSON file.
* `credentials_json`: Path to the credentials JSON file.

**Methods**:

* `run()`: Runs the inspection checks on the specified model.
* `_get_checks(inspector)`: Retrieves the list of checks to perform on the model.

**Example**:

```python
inspector = InspectService(dir="/path/to/model", model="eosxxxx")
results = inspector.run()
```

***

### Class: `CheckService`

Service for performing various checks on the model.

**Parameters:**

* `logger`: Logger for logging messages.
* `model_id`: Identifier of the model.
* `dest_dir`: Destination directory for storing model-related files.
* `dir`: Directory where the model repository is located.
* `ios`: Instance of `IOService` for handling input/output operations.

**Methods**:

* `check_files()`: Checks the existence of required files for the model.
* `check_information(output)`: Performs various checks on the model information.
* `check_single_input(output, run_model, run_example)`: Checks if the model can run with a single input without error.
* `check_example_input(output, run_model, run_example)`: Checks if the model can run with example input without error.
* `check_consistent_output(run_example, run_model)`: Checks if the model produces consistent output.

**Example**:

```python
check_service = CheckService(logger=logger, model_id="eosxxxx", dest_dir="/path/to/dest", dir="/path/to/dir", ios=ios)
check_service.check_files()
```

***

### Class: `IOService`

Service for handling input/output operations related to model testing.

**Parameters:**

* `logger`: Logger for logging messages.
* `dest_dir`: Destination directory for storing model-related files.
* `model_path`: Path to the model.
* `bundle_path`: Path to the model bundle.
* `bentoml_path`: Path to the BentoML files.
* `model_id`: Identifier of the model.
* `dir`: Directory where the model repository is located.

**Methods**:

* `_run_check(check_function, data, check_name, additional_info)`: Runs a check function and logs the result.
* `_generate_table(title, headers, rows, large_table, merge)`: Generates a table with the given title, headers, and rows.
* `get_model_type(model_id, repo_path)`: Gets the type of the model based on the repository contents.
* `get_file_requirements()`: Gets the list of required files for the model.
* `read_information()`: Reads the information file for the model.
* `print_output(result, output)`: Prints the output of a result.
* `get_conda_env_size()`: Gets the size of the Conda environment for the model.
* `calculate_directory_size(path)`: Calculates the size of a directory.
* `get_directories_sizes()`: Gets the sizes of the model directory and the Conda environment directory.

**Example**:

```python
ios = IOService(logger=logger, dest_dir="/path/to/dest", model_path="/path/to/model", bundle_path="/path/to/bundle", bentoml_path="/path/to/bentoml", model_id="eosxxxx", dir="/path/to/dir")
ios.read_information()
```
