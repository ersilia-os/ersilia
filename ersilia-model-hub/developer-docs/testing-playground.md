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
The Testing Playground only works with Linux systems, and is oriented towards Ersilia developers with experience in our tools.
{% endhint %}

## TL:DR

To use the Test Playground, you need to have Ersilia installed in test mode and the Ersilia repository cloned into your local.

```bash
conda create -n ersilia python=3.12
conda activate ersilia
git clone https://github.com/ersilia-os/ersilia.git
pip install -e .[test]
```

The Playground runs on a `nox` environment, isolated from the local source. Each time we will activate a `nox` session from the `test/playground` folder in Ersilia:

<pre class="language-bash"><code class="lang-bash">cd ersilia/test/playground
<strong>nox -s execute # to start the session and run commands
</strong>nox -s clean #to clean up all folders and session files created by nox
</code></pre>

The `nox` command `nox -s execute` will initiate a session. A few built-in flags can be passed to it:

<table><thead><tr><th width="117">Flag</th><th width="170">Status</th><th>Description</th></tr></thead><tbody><tr><td>-s</td><td>Required</td><td>Specifies the session to be run by Nox.</td></tr><tr><td>-p</td><td>Not required</td><td>Specifies the python environment. If none is specified, it will test everything in py3.8 to py3.12. More than one environment can be specified simply using <code>nox -s execute -p 3.8 3.9</code>.</td></tr><tr><td>-fb</td><td>Not required</td><td>Used to change python backends (conda, mamba, micromamba, virtualenv, venv, uv, none). Defaults to conda.</td></tr><tr><td>-v</td><td>Not required</td><td>Verbose output printed in the terminal.</td></tr></tbody></table>

{% hint style="warning" %}
The first time you use Nox in your system you will be required to grant sudo privileges so that actions like DockerHub activation can be performed
{% endhint %}

{% hint style="info" %}
Please note that nox separates built-in arguments and custom flags using `--`. Hence, all arguments not referred in this table will need to follow the following structure:

`nox [built-in flags] -- [custom flags]`  &#x20;
{% endhint %}

### Structure overview

The idea behind the playground is to cover all sorts of tests we might want to do on any model, more extensively than the test command itself. Therefore, it is a highly customizable functionality and only addressed to Ersilia developers.

Nox will create an isolated environment and store the files used for testing under `~/eos/playground/files` and the logs generated under `~/eos/playground/logs`. Those will be eliminated with the `nox -s clean` command.&#x20;

The playground is adapted to many use cases, for example:

* Test several models on python 3.12 fetching them from github
* Test one single model across all python versions fetching from dockerhub
* Evaluate if a model seems to have gotten slower
* ...

Below we describe the flags you can use to combine all these custom-made tests.

## Ersilia playground flags

If no flag is specified, the command `nox -s execute` will run its default test: it will try to fetch, serve and run the model `eos3b5e` (molecular weight) in all python environments from 3.8 to the latest maintained by Ersilia. The model will be fetched from github by default and a shallow test will also be performed (see [model test](model-tester.md) section). More granularity can be specified using the built-in flags (for python versions and environments) as well as the flags specific to the Ersilia CLI:

<table><thead><tr><th width="133">Flag</th><th width="147">Default</th><th width="238">Description</th><th>Example</th></tr></thead><tbody><tr><td>--cli</td><td>all</td><td>Specifies ersilia commands to run in order (fetch, serve, run, catalog, example, test, close, delete). Default is all, which executes commands in this order: "fetch", "serve", "run", "close", "catalog", "example", "delete", "test".</td><td><code>nox -s execute -- --cli fetch serve run</code><br><code>nox -s execute -- --cli run</code> </td></tr><tr><td>--fetch</td><td>--from_github</td><td>Fetches models from sources (from_github, from_dockerhub, from_s3, version)</td><td><p><code>nox -s execute -- --fetch from_s3</code></p><p><code>nox -s execute -- --fetch from_dockerhub version dev</code></p></td></tr><tr><td>--run</td><td>None</td><td>Run a model. It will try to fetch it from_dockerhub as this is Ersilia's default (depends on the activate_docker flag, see below). Best combined with the input and output flags (see below)</td><td>nox <code>-s execute -- --run</code></td></tr><tr><td>--example</td><td>["-n", 10, "--random"]</td><td>Generates example input for a model (-n, --random, -f). If we specify a </td><td><code>nox -s execute -- --example -n 10 random/predefined -c -f example.csv</code></td></tr><tr><td>--catalog</td><td>["--more", "--local", "--as-json"]</td><td>Retrieves model catalog from local or hub.</td><td><code>nox -s execute -- --catalog hub</code></td></tr><tr><td>--test</td><td>["--shallow", "--from_github"]</td><td>Tests models at different levels (shallow and deep) and from_github, from_dockerhub or from_s3.</td><td><code>nox -s execute -- --test deep from_dockerhub/from_s3/from_github</code></td></tr><tr><td>--delete</td><td>None</td><td>Used to delete models. It has only one flag: all</td><td><p></p><p><code>nox -s execute -- --delete all</code><br></p></td></tr></tbody></table>

### **Additional flags for Ersilia's CLI**

<table><thead><tr><th width="133">Flag</th><th width="147">Default</th><th width="238">Description</th><th>Example</th></tr></thead><tbody><tr><td>--outputs</td><td>[results.{csv, json, h5}]</td><td>This is used with run command and used to specify output file types. Note that we only specified the file name, the path will be automatically set to ~/eos/playground/files/{file_name.{csv, json, h5}}</td><td><code>nox -s execute -- --outputs result.csv result.h5</code></td></tr><tr><td>--input_types</td><td>List of (str, list, csv)</td><td>This is also used with run command to define input formats (str, list, csv).</td><td><code>nox -s execute -- --input_types str list csv</code></td></tr><tr><td>--runner</td><td>single</td><td>Specifies execution mode (single, multiple). The single mode is used to execute commands using one model whereas the multiple mode will use multiple models to execute the given commands.</td><td><code>nox -s execute -- --runner multiple</code></td></tr><tr><td>--single</td><td>eos3b5e</td><td>Used to specify or override the default model ID used for single running mode.</td><td><code>nox -s execute -- --single model_id</code></td></tr><tr><td>--multiple</td><td>[eos5axz, eos4e40, eos2r5a, eos4zfy, eos8fma]</td><td>Used to specify or override the default model IDs used for multiple running mode.</td><td><code>nox -s execute -- --multiple model_id1 model_id2</code></td></tr></tbody></table>

### General setting flags

<table><thead><tr><th width="171">Flag</th><th width="147">Default</th><th width="238">Description</th><th>Example</th></tr></thead><tbody><tr><td>--activate_docker</td><td>true</td><td>Activates or deactivates Docker. It allows to test for example if autofetcher will decide not to fetch from Docker (default) when Docker is not active</td><td><code>nox -s execute -- --activate_docker false</code></td></tr><tr><td>--log_error</td><td>true</td><td>Enables or disables logging of errors as file, which will be stored in ~/eos/playground/logs/. Each command failures will create a standalone file, with datetime on it in a string format. For instance catalog_20250129_145802.txt</td><td><code>nox -s execute -- --log_error false</code></td></tr><tr><td>--silent</td><td>true</td><td>Enable or disable logs from ersilia command execution</td><td><code>nox -s execute -- --silent false</code></td></tr><tr><td>--show_remark</td><td>false</td><td>Displays a remark column in the final execution summary table that allows to quickly see if the ersilia commands are executed successfully</td><td><code>nox -s execute -- --show_remark true</code></td></tr><tr><td>--max_runtime_minutes</td><td>10</td><td>Sets the maximum execution time for a run command, to test model speed if seems to be slow.</td><td><code>nox -s execute -- --max_runtime_minutes 5</code></td></tr><tr><td>--num_samples</td><td>10</td><td>Sets the sample size to create input for <code>run</code> command.</td><td><code>nox -s execute -- --num_samples 5</code></td></tr></tbody></table>
