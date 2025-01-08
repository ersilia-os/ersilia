---
description: >-
  The Testing Playground is a new testing system added to the Ersilia testing
  pipeline. It aims to advance the testing pipeline using the ersilia CLI
  commands. It provides a flexible and robust fram
layout:
  title:
    visible: true
  description:
    visible: false
  tableOfContents:
    visible: true
  outline:
    visible: true
  pagination:
    visible: true
---

# Testing Playground

The **Testing Playground** provides a flexible and robust framework for validating and profiling CLI commands that are being used for managing Ersilia models from various sources, such as GitHub, DockerHub, and local directories.

{% hint style="warning" %}
The Testing Playground only works with Linux systems.
{% endhint %}

## Usage

To use the Test Playground, you need to have Ersilia installed in test mode and the Ersilia repository cloned into your local, as the files will be saved in the Ersilia directory

```bash
conda create -n ersilia python=3.11
conda activate ersilia
pip install ersilia[test]
cd ersilia
nox -f test/playground/noxfile.py -s setup
```

### **Tests available**

1. From GitHub: the model is fetched from its source directory and installed in an independent conda environment --> `nox -f test/playground/noxfile.py -s test_from_github`
2. From DockerHub: the model is fetched via its Docker Image and run inside a container -->
3. AutoFetcher: the test will automatically decide which fetching system to use
4. Fetch multiple models:
5. Serve multiple models:
6. Conventional run

### **Playground Key Features**

1. **Enhanced Testing for CLI Commands**:
   * Supports five CLI commands currently: `delete`, `fetch`, `serve`, `run`, `close`.
   * Ensures accurate validation of model and CLI behavior and performance.
2. **Customizable Rules**:
   * Users can define custom rules using the `rules.py` file, leveraging the `Registry Design Pattern` to make sure single instance of rule instantiated.
   * Rules can specify expected behavior and perform checks after each command execution.
3. **Resource Profiling**:
   * After a single check performed, its memory usage, execution time, and check pass or fail of our expected behavior will be captured.
   * Logs errors and exceptions to assist with debugging.
4. **Pytest Integration**: -Every checks are using pytest and expecetd to pass several assertions.
   * Collects test results in a shared list (`shared.py`) for aggregation since pytest does not allow a direct or nice reporting display in the pytest code.
   * A global config file defined in `conftest.py` at the root of project directory which displays results in a **Rich** table with memory, runtime, and status information with pytest printhook (`pytest_terminal_summary(terminalreporter, exitstatus, config)`).
5. **Configurable Execution**:
   * The playground has `config.yml` file that supports several configureable options including single or multiple commdand execution of models, selecting command type, and more.
6. **Nox Integration**:
   * The playground is powered by `nox`. `nox` run tests in isolated or existing Python environments using `nox` sessions. All the sessions are defined in the `noxfile.py`.

### **How It Works?**

1.  **PLayground folder structure**

    The playground folder structure are given below, one folder that wont't be explained anyhere in detail is \`files\`. It is simply used for input, output, and error logs files to be stored, to make the playgroud folder a bit cleaner.\


    ```plaintext
    playground
    ├── commands.py
    ├── config.yml
    ├── files
    │   ├── error_log_20241122_133958.txt
    │   ├── input.csv
    │   ├── result.csv
    │   └── sample_error_log.txt
    ├── __init__.py
    ├── noxfile.py
    ├── rules.py
    ├── shared.py
    └── utils.py
    ```
2.  **Defining Rules**:

    The user creates rules in \`rules.py\` to validate some expected behavior. These rules are executed after each CLI command is executed. For instance we want to check folder eixtence after we fetch a model, this rule is going to get executed after the fetch command is done. Our expected behaviour is \`\`\`True\`\`\` in this case, we want the folder of the model exists. The rule structure mainly contains \`rule registry\`, \`the rule it self\` and a function to access registered rules (\`get\_rules\`).

    Here is rule registry decorator looks like:



    ```python
     RULE_REGISTRY = {}
    ```

    ```python
     def register_rule(name):
         def decorator(cls):
             RULE_REGISTRY[name] = cls
             return cls

         return decorator
    ```

    \


    Here is an example rule for checking folder existence for the fetched model. When the rule is defined the minimum structure has to be satified specially the exception or assertion handling part.\


    ```python
     @register_rule("folder_exists")
     class FolderExistsRule(CommandRule):
         def __init__(self):
             pass

         def check(self, folder_path, expected_status):
             actual_status = Path(folder_path).exists() and any(
                 Path(folder_path).iterdir()
             )
             if actual_status != expected_status:
                 raise AssertionError(
                     f"Expectation failed for FolderExistsRule: "
                     f"Expected folder to {'exist' if expected_status else 'not exist'}, "
                     f"but it {'exists' if actual_status else 'does not exist'}."
                 )
             return {
                 "name": f"Folder exists at {folder_path}",
                 "status": actual_status,
             }
    ```


3.  **Run Commands**:

    The rule defined above will be registed in the \`get\_rule\` function below. The design pattern for the rule definition ensures that onc instance of a single rule will be created and accessed.\


    ```python
     def get_rule(rule_name, *args, **kwargs):
         rule_class = RULE_REGISTRY.get(rule_name)
         if not rule_class:
             raise ValueError(f"Rule '{rule_name}' is not registered.")
         return rule_class().check(*args, **kwargs)
    ```



    In the `commands.py` we import our rules and put them in the `apply_rules` function after the fetch section for our case and lets assume that we fetched from `--from_dockerhub`:\


    ```python
     def apply_rules(command, description, dest_path, repo_path, config):
         checkups = []
         try:
             if description == "fetch":
                 if from_dockerhub:
                     checkups.append(
                         get_rule(
                             rule_name="folder_exists",
                             folder_path=dest_path,
                             expected_status=True,
                         )
                     )
         except Exception as rule_error:
             handle_error_logging(
                 command,
                 description,
                 rule_error,
                 config,
                 checkups
             )
             pytest.fail(f"Rule exception occurred: {rule_error}")

         return checkups
    ```



    The \`apply\_rules\` function will then be executed in \`execute\_commands\` function like the snippet below. The \`create\_compund\_input\_csv\` gets executed at the start of command of execution, helps to generate input example at run time to enable auto input example file generation at github workflow. But locally we can use our own file by specifying at the \`config.yml\` with \`input\_file\` section and just remove this function. We may need some option for that in next update:\


    ```python
         def execute_command(command, description=None, dest_path=None, repo_path=None):
             create_compound_input_csv(config.get("input_file"))

             start_time, max_memory, success, result, checkups, = time.time(), 0, False, "", [] 

             proc = psutil.Popen(
                 command, 
                 stdout=subprocess.PIPE, 
                 stderr=subprocess.PIPE
             )

             try:
                 while proc.poll() is None:
                     max_memory = max(max_memory, proc.memory_info().rss / (1024 * 1024))
                     time.sleep(0.1)

                 success = proc.returncode == 0
                 stdout, stderr = proc.communicate()

                 if success:
                     result = stdout.decode()
                 else:
                     result = stderr.decode()

             except Exception as e:
                 proc.kill()
                 result = str(e)

                 if config.get("log_error", False):
                     handle_error_logging(
                         command,
                         description,
                         result,
                         config
                     )

                 pytest.fail(
                     f"{description} '{' '.join(command)}' failed with error: {result}"
                 )

             checkups = apply_rules(
                 command, 
                 description, 
                 dest_path, 
                 repo_path,
                 config
             )
    ```



    Then after this, we have a pytest function that does the tests for the commands. Here is the pytest code snippet for that.



    ```python
     @pytest.mark.parametrize("command_name", get_command_names(model_ids[0], cli_type, config))
     def test_command(model_id, command_name):
         commands = get_commands(model_id, config)
         command = commands[command_name]
         dest_path = base_path / "dest" / model_id
         repo_path = base_path / "repository" / model_id

         success, time_taken = execute_command(
             command,
             description=command_name,
             dest_path=dest_path,
             repo_path=repo_path,
         )

         assert success, f"Command '{command_name}' failed for model ID {model_id}"
    ```



    Note that the \`get\_command\_names\` function provides all the five commands and structured and appropriate way to efficiently run them on pytest test function in parameterized way. This helper function is defined in \`utils.py\` file.\

4.  **Capture Results and Display**:

    To be able to share data from a \`commands.py\` to a pytest print hook defined at the root in \`conftest.py\`, the \`shared.py\` was created. This file simple defined an empty list \`result = \[]\`, to store every necessary info from a pytest execution as below:\


    ```python
     from .shared import results
     results.append(
         {
             "command": " ".join(command),
             "description": description,
             "time_taken": f"{(time.time() - start_time) / 60:.2f} min",
             "max_memory": f"{max_memory:.2f} MB",
             "status": status_text,
             "checkups": checkups,
             "activate_docker": docker_activated,
             "runner": config.get("runner"),
             "cli_type": config.get("cli_type"),
         }
     )
    ```

    \


    Then aggregated test outcomes in a `shared.py` list will be presented in a Rich table within `conftest.py` then to the terminal.



    Then the final thing is to call and run the pytest command file \`commands.py\` in the nox session runner \`noxfile.py\`. So by defualt ersilia has multiple selected models(defined in \`config.yml\` at \`model\_ids\` section) and we run several tests on them. Those tests include, for a runner type \`single\` are \`auto fetcher\`: to fetch models automatically from a recommened source, \`fetch from github\`, \`fetch from docker hub\`, \`serve the fetched model\` and run either \*\*Conventional Runner\*\* (its a runner that does a detail checkups on the input files schema but slower) and \*\*Standard Runner\*\* (which accepts standard input file and does a few checkups but faster) will be run. For runner type \`multiple\` fetching from a dockerhub for all \`model\_ids\` at once and serve them in batch. Before that an example session code given for \`fetching from github\` and runner type \`single\` as below:\


    ```python
     @nox.session(venv_backend="conda", python=get_python_version())
     def setup(session):
         # to install dependency
         install_dependencies(session)
         # to install ersilia either from a source or from dev environment
         setup_ersilia(session)
     
     # Sample nox session to execute fetch fro github
     @nox.session(venv_backend="conda", python=get_python_version())
     def test_from_github(session):
         install_dependencies(session)
         logger.info(
             f'CLI test for model: {config.get("model_id")} and {config.get("fetch_flags")}'
         )
         session.run("pytest", "commands.py", "-v", silent=False)
    ```



    If we need to create another session we can update the config file using `update_yaml_values` such as changing the runner type or other things. To run a nox session, we need to go to the playground folder and type the following. To only run a custom session we type:



    ```bash
    nox -s setup test_from_github

    ```



    To run all sessions

    ```bash
    nox
    ```



    Or we can run them from any folder but by specifying the where the nox file existed using `-f /path/to/noxfile.py`. Below is an example of minimal terminal output



    ```plaintext
     ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
     ┃ Command                       │ Status │ Runtime    │ Memory Usage   ┃
     ┣━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ ┫
     ┃ ersilia fetch eos3b5e         │ PASSED │ 1.23 min   │ 128.00 MB      ┃
     ┃ ersilia serve eos3b5e         │ PASSED │ 0.45 min   │ 64.00 MB       ┃
     ┃ ersilia run -i input.csv      │ FAILED │ 0.78 min   │ 96.00 MB       ┃
     ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
     
    ```

### **Configuration**

The `config.yml` file specifies various options for customizing the playground structured as below:

```yaml
activate_docker: true           # To disable and enable docker at test time
cli_type: all                   # Allow us to specify the commnad we want to run
delete_model: true              # Allow us to delete fetched model if previously existed
fetch_flags: --from_github      # To fetch from a specified source
input_file: files/input.csv     # Input file path
log_error: true                 # Enables logging the error as a file
log_path: files/error_log       # Error logging path, better to in files folder
max_runtime_minutes: 10         # Expected runtime for run command
model_id: eos3b5e               # single model id
model_ids:                      # Defines multiple mode ids
  - eos3b5e                     |
  - eos4e40                     |    
  - eos9gg2                     |
output_file: files/result.csv   # Path for output file
output_redirection: false       # Sometimes we need our output to be redirected to json file, this enables that
overwrite_ersilia_repo: false   # This allow us to over our existed isolated venv
python_version: 3.10.10         # ---
run_flags: ''                   # flags that can be add in the run command like '-o iles/result.csv'
runner: single                  # enables single or multiple commands to be run at once
serve_flags: ''                 # flags at serve time like "--track_runners"
use_existing_env: true          # enables to use our existing venv(base conda venv for now)
```

### **Github action workflow intergation**

In addition to using the playground locally, its also a part of in Ersilia's CI/CD workflow. Its integrated at \`test\_and\_cleanup.yml\` workflow file. For a single runner, the sessions are independent so that they were implemented to be parallely executed. But for the multiple runner they are dependent on each other (we can not serve before fetch gets completed), they were implemented to run one after the other. Example codes below:

```yaml
strategy:
    matrix:
    session:
        - setup
        - test_from_github
        - test_from_dockerhub
        - test_auto_fetcher_decider
        - test_conventional_run

- name: Run CLI Test Default
    run: |
    source activate            
    nox -f test/playground/noxfile.py -s ${{ matrix.session }}
```
