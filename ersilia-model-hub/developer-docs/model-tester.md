# Model Tester

### Overview

This functionality provides tools and tests for running and validating model operations that are underdevelopment locally or even deployed models remotely. Models in ersilia are packed with two methods, `bentoml` and `ersilia-pack`. This CLI supports models packed with both types of mehtods. Checks includes:

* **High-Level Checks:**
  * Model Size Calculation: Checks the size of the `venv` of the model and othe size for the directory itself. In the first case it helps to note the package sizes necessary to run the model and the second case weighs the files (this includes model checkpoints) necessary to run the model.
  * Model Tasks: the model supports a variety of tasks, including classification, regression, generative modeling, representation learning, similarity assessment, clustering, and dimensionality reduction, providing flexibility to address diverse machine learning and data analysis needs.
  * Input and Output Shapes: for input shape to ensure it matches one of the supported formats: a single entity, a pair of entities, a list of entities, a pair of lists, or a list of lists. These checks are performed at runtime to guarantee compatibility and prevent errors, ensuring smooth and efficient processing of data. On the other hand the model's output will conform to one of the supported shapes: a single entity, a list of entities, a flexible list (with varying sizes), a matrix, or a serializable object. This ensures that the output format aligns with the expected structure for downstream processing and facilitates seamless integration.
  * Model Outputs Types: The model produces outputs in various formats, including boolean values, compounds, descriptors, distances, experimental values, images, other values, probabilities, proteins, scores, and text, offering versatility to accommodate a wide range of applications.
  * File Integrity Checks: this checks a models packed with both methods have a files necessary to be there and flags if missed.
* **Detailed Inspection with `--level/-l deep` Flag:**
  * A deeper inspection can now be performed to include:
    * Testing with a single input.
    * **Model consistency output:** Comparing outputs from running `run.sh` and the Ersilia `run` command. So the run.sh create an isolated vertual environment and generate results from the shell. The other one uses `ersilia` CLI to fetch, serve and run models. The output from both method will be compared using `check_consistent_output` function. Checks the string and numerical output consistency will be calculated. The scores include `rmse`, `spearmanr` (Spearman correlation coefficient with associated p-value) etc...
    * Computational performance assessment: After serving the models, they get executed for 1, 50, and 100 input and performance using `wall clock` will be recorded and reported.
    * Complete Metadata: this check ensures the metadata file's completeness by verifying required fields and validating URLs. The required fields are `["Publication", "Source Code", "S3", "DockerHub"]`. The function first constructs the metadata file URL based on the package type (`bentoml` or otherwise) and checks if the file exists. If the file is missing or cannot be fetched/parsed, it returns a failure result. It then identifies missing fields and invalid URLs within the metadata, appending relevant details. Additionally, the function attempts to parse the metadata file using `RepoMetadataFile.read_information`, catching and reporting any exceptions. If all checks pass, it confirms the metadata is complete; otherwise, it consolidates the details of the issues and returns a failure result.
    * Repository structure validation, detection of extra/unnecessary files, and more.

#### Mechanism

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

Command line example

```bash
$ ersilia test model_id -d /path/to/model --remote --inspect --remove -l deep
```
