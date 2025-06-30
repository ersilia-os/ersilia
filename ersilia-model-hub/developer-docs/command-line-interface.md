---
description: Command Line Interface to interact with the Ersilia Model Hub
---

# Command line interface

### 1. Fetch Command

Fetch an Ersilia model to run it locally. Ersilia defaults to using Docker images of the models to avoid system-related problems, your Docker engine needs to be active to allow this. Docker images run containers completely isolated from your system so they will not interfere with any installation you have. Otherwise, models can be fetched from source, i.e, from a local, cloned or downloaded repository. In this case, a conda environment is created to install the necessary packages for the model to run.  When a model is re-fetched, Ersilia deletes all artifacts related to it, but for the sake of efficiency, the user is prompted to override this behavior.

<table><thead><tr><th width="256">Option / Argument</th><th width="102">Type</th><th width="106">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>&#x3C;MODEL_ID></code> (argument)</td><td>String</td><td><em>Required</em></td><td>The model ID to be fetched. Cn either be the eos identifier or the slug identifier.</td></tr><tr><td><code>--from_dir</code></td><td>String</td><td><code>None</code></td><td>Fetch the given model, from source, using a local path provided on the user's system. The path can be relative or absolute. </td></tr><tr><td><code>--from_github</code></td><td>Flag</td><td><code>False</code></td><td>Fetch the given model, from source, by cloning the GitHub repository of the model on the user's system.</td></tr><tr><td><code>--from_s3</code></td><td>Flag</td><td><code>False</code></td><td>Fetch the given model, from source, using the model's code archive maintained on Ersilia's S3 storage servers.</td></tr><tr><td><code>--from_dockerhub</code></td><td>Flag</td><td><code>True</code></td><td>Fetch the given model using its Docker image maintained in Ersilia's public DockerHub registry.</td></tr><tr><td><code>--version</code></td><td>String</td><td><code>latest</code></td><td>Version of the model to fetch, used when fetching from DockerHub. Only for testing purposes.</td></tr><tr><td><code>--from_hosted</code></td><td>Flag</td><td><code>False</code></td><td>Fetch the given model, using the URL where the model is hosted. This only creates a basic folder structure for the model, and the model is not actually downloaded.</td></tr><tr><td><code>--with_fastapi</code></td><td>Flag</td><td><code>True</code></td><td>This is a developer specific command and most users will not need to use it. This command fetches a model from source, and packs it using Ersilia Pack, using FastAPI.  During model contribution and troubleshooting, it is often useful to fetch a model again, however one may not want to rebuild the model's environment, in which case it can be useful to combine this flag with the <code>--reuse</code> flag. </td></tr><tr><td><code>--overwrite/--reuse</code></td><td>Boolean</td><td><code>True</code> (overwrite)</td><td>This flag is useful from an efficiency perspective. When fetching a model from source, a conda environment is created, which can be retained and reused, or overwritten, when re-fetching a model.</td></tr></tbody></table>

**Example:**

```bash
$ ersilia fetch eosxxxx --from_github/--from_dockerhub/--from_s3/--from_hosted
```

```
$ ersilia fetch eosxxxx --from_github/--from_dockerhub/--from_s3/--from_hosted --overwrite
```

***

### 2. Serve Command

This command serves a specified model as an API.

<table><thead><tr><th width="237">Option / Argument</th><th width="96">Type</th><th width="113">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>&#x3C;MODEL_ID></code> (argument)</td><td>String</td><td><em>Required</em></td><td>The model ID to be served. Cn either be the eos identifier or the slug identifier.</td></tr><tr><td><code>--port</code> / <code>-p</code></td><td>Integer</td><td><code>None</code></td><td>The port to use when creating a model server. If unspecified, Ersilia looks for empty ports to use on the user's system.</td></tr><tr><td><code>--track</code> / <code>-t</code></td><td>Flag</td><td><code>False</code></td><td>Whether the model's runs should be tracked to monitor model and system performance. This telemetry data can be uploaded to S3 if the appropriate credentials are provided. Currently only Ersilia developers and other Ersilia tools have that privilege.</td></tr><tr><td><code>--tracking-use-case</code></td><td>String</td><td><code>local</code></td><td>If <code>--track</code> is true, this command allows specification of the tracking use case. Current options are: <code>local</code>, <code>hosted</code>, <code>self-service</code> and <code>test</code>.</td></tr><tr><td><code>--enable-local-cache/--disable-local-cache</code></td><td>Flag</td><td><code>True</code></td><td>Toggle Redis-based local caching on or off. If it enabled the results from model APIs will be cached for 7 days.</td></tr><tr><td><code>--max-cache-memory-frac</code></td><td>Float</td><td><code>None</code></td><td>Sets the maximum fraction of memory to use by Redis for caching. Recommened value <code>0.2-0.7.</code></td></tr><tr><td><code>--local-cache-only</code></td><td>Bool</td><td><code>False</code></td><td>Specifies to fetch stored model results from local cache. The local caching system is powered by Redis.</td></tr><tr><td><code>--cloud-cache-only</code></td><td>Bool</td><td><code>False</code></td><td>Specifies to fetch stored model results from cloud cache. This allows to fetch model precalculated results in csv file in Ersilia model output format.</td></tr><tr><td><code>--cache-only</code></td><td>Bool</td><td><code>False</code></td><td>Specifies to fetch stored model results from both local and cloud cache. More details are give in a dump CLI.</td></tr></tbody></table>

**Examples:**

```bash
$ ersilia serve eosxxxx --no-cache --port 12450
```

```
$ ersilia serve eosxxxx --max-cache-memory-frac 0.5
```

***

### 3. Run Command

This command runs a specified model with given inputs and outputs in the shell where the model gets served.

<table><thead><tr><th width="189">Option / Argument</th><th width="97">Type</th><th width="118">Default Value</th><th width="272">Description</th><th data-type="checkbox">Required</th></tr></thead><tbody><tr><td><code>-i, --input</code></td><td>String</td><td><code>None</code></td><td>Pass the input for the model as a .csv file with only one column with header.</td><td>true</td></tr><tr><td><code>-o, --output</code></td><td>String</td><td><code>None</code></td><td>Specify the file to save the model predictions. This can be a CSV or H5.</td><td>true</td></tr><tr><td><code>-b, --batch_size</code></td><td>Integer</td><td><code>100</code></td><td>Specify the batch size for generating model predictions. By default, Ersilia works with batch size of 100 inputs.</td><td>false</td></tr><tr><td><code>--as_table</code> / <code>-t</code></td><td>Flag</td><td><code>False</code></td><td>Print the model predictions in an ASCII table on the terminal.</td><td>false</td></tr></tbody></table>

***

**Example:**

```bash
$ ersilia run -i input.csv -o output.csv --batch_size 1000 --as_table
```

### 4. Example Command

This command can sample inputs for a given model. By default, the command expects a model to be served in the current shell, otherwise a model needs to be explicitly specified. By default, five examples will be sampled for the model and displayed on the user's terminal. This behavior can be customized with relevant flags.

<table><thead><tr><th width="242">Option / Argument</th><th width="110">Type</th><th width="137">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><code>None</code></td><td>Optional model ID for which to generate examples, if no model is already served.</td></tr><tr><td><code>--n_samples</code>, <code>-n</code></td><td>Integer</td><td><code>5</code></td><td>Specify the number of example inputs to generate for the given model.</td></tr><tr><td><code>--file_name</code>, <code>-f</code></td><td>String</td><td><code>None</code></td><td>File name where the examples should be saved.</td></tr><tr><td><code>--simple</code> / <code>--complete</code></td><td>Flag</td><td><code>--simple</code> </td><td>Simple inputs only contain the SMILES, while complete inputs also include InChIKey and the molecule's name.</td></tr><tr><td><code>--random</code> / <code>--predefined</code></td><td>Flag</td><td><code>--random</code></td><td>If the model source contains an example input file, when the predefined flag is set, then inputs are sampled from that file. Only the number of samples present in the file are returned, especially if --n_samples is greater than that number. By default, Ersilia samples inputs randomly.</td></tr><tr><td><code>--deterministic</code></td><td>Flag</td><td><code>False</code></td><td>Used to generate examples data deterministically instead of random sampling. This allows when every time you run with example command with this flag you get the same types of examples. </td></tr></tbody></table>

**Example:**

```bash
$ ersilia example -n 10 --file_name examples.csv
```

***

### 5. Catalog Command

List the catalog of Ersilia models. By default, this command only shows the catalog of models available locally on a user's system, in an ASCII table on the user's terminal.

<table><thead><tr><th width="244">Option / Argument</th><th width="78">Type</th><th width="135">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>--hub/--local</code></td><td>Flag</td><td><code>--local</code></td><td>Toggle between showing models available in the model hub or locally.</td></tr><tr><td><code>--file_name</code>, <code>-f</code></td><td>String</td><td><code>None</code></td><td>Write the requested catalog to a given file.</td></tr><tr><td><code>--more</code> / <code>--less</code></td><td>Flag</td><td><code>--less</code></td><td>Print more or less information about the catalog. With less (default) only the <em>identifier</em> and <em>slug</em> are displayed. Otherwise, the following fields are displayed: <em>identifier, slug, title, task, output, output shape,</em> and <em>model source</em>.  These flags work with both local and hub context.</td></tr><tr><td><code>--card</code></td><td>Flag</td><td><code>&#x3C;MODEL_ID></code></td><td>Displays the card for the specified model identifier. The requested model does not need to be available locally. </td></tr><tr><td><code>--as-json /--as-table</code></td><td>Flag</td><td><code>False</code></td><td>Displays the catalog in JSON format on the user's terminal. This may be helpful for downstream processing task. The default behavior to display the catalog in a tabular format. These flags can be used both in the local and hub context.</td></tr></tbody></table>

**Examples:**

```bash
$ ersilia catalog
```

```
$ ersilia catalog --more --hub
```

```
$ ersilia catalog --card eosxxxx --as-json
```

### 6. Delete Command

This command deletes a specified model from the local storage. You can delete a specific model by providing its ID or remove all locally available models using the `--all` flag.

<table><thead><tr><th width="174">Option / Argument</th><th width="74">Type</th><th width="118">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><em>Optional</em></td><td>The model ID to delete.</td></tr><tr><td><code>--all</code></td><td>Flag</td><td><code>False</code></td><td>When provided, deletes all locally available models.</td></tr></tbody></table>

**Examples:**

```bash
$ ersilia delete eosxxxx
```

```bash
$ ersilia delete --all
```

***

### 7. Close Command

This command closes the current session of the served model and cleans up any associated resources.

| Option / Argument         | Type | Default Value | Description                                |
| ------------------------- | ---- | ------------- | ------------------------------------------ |
| _(No additional options)_ |      |               | Just run the command to close the session. |

**Example:**

```bash
$ ersilia close
```

***

### 8. Dump Command

This command used to fetch cached precalculation from the cloud and the local Redis system depending on the flag specified in the serve command [#id-2.-serve-command](command-line-interface.md#id-2.-serve-command "mention"). <mark style="color:red;">Note that fetching sample size of <= 50,000 from a cloud is not recommended as it is expensive.</mark>

<table><thead><tr><th width="196">Option / Argument</th><th width="92">Type</th><th width="78">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>-n/--n_samples</code></td><td>Integer</td><td>-1</td><td>The sample size to for the cache we want fetch either from a cloud or from a local. Its is <code>-1</code> by default which indicates fetch all stored cached results.</td></tr><tr><td><code>-o/--output/output</code></td><td>String</td><td>None</td><td>The path of the csv file to store the fetched cached results. Currently this CLI only supports csv as its output.</td></tr></tbody></table>

**Example:**

```bash
$ ersilia dump -o my_output.csv
```

```
 ersilia dump -n 100000 -o my_output.csv
```

### 9. Test Command

This command tests a model and produces performance metrics, allowing you to specify various options such as the source from which to fetch the model, and the depth of the testing process.

<table><thead><tr><th width="182">Option / Argument</th><th width="80">Type</th><th width="109">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>&#x3C;MODEL_ID></code> (argument)</td><td>String</td><td><em>Required</em></td><td>The model ID to test, either eos identifier or slug.</td></tr><tr><td><code>--from_dir</code></td><td>String</td><td><code>None</code></td><td>Local directory where the model is stored.</td></tr><tr><td><code>--from_github</code></td><td>Flag</td><td><code>True</code></td><td>Fetch the model directly from GitHub.</td></tr><tr><td><code>--from_dockerhub</code></td><td>Flag</td><td><code>False</code></td><td>Force fetching the model from DockerHub.</td></tr><tr><td><code>--from_s3</code></td><td>Flag</td><td><code>False</code></td><td>Force fetching the model from an AWS S3 bucket.</td></tr><tr><td><code>--version</code></td><td>String</td><td><code>latest</code></td><td>Specify the version of the model to fetch, especially useful when fetching from DockerHub.</td></tr><tr><td><code>--inspect</code></td><td>Flag</td><td><code>False</code></td><td>Run inspection tests for further model diagnostics.</td></tr><tr><td><code>--surface</code></td><td>Flag</td><td><code>False</code></td><td>Additional flag for testing, similar to deep checks.</td></tr><tr><td><code>--shallow</code></td><td>Flag</td><td><code>False</code></td><td>Run shallow tests (e.g., container size, output consistency).</td></tr><tr><td><code>--deep</code></td><td>Flag</td><td><code>False</code></td><td>Run deep tests (e.g., computational performance checks).</td></tr><tr><td><code>--report_path</code></td><td>String</td><td><code>&#x3C;MODEL_ID>_test.json</code></td><td>Specify the file path to output the test report in JSON format. Automatically set to the root from where the command is called.</td></tr><tr><td><code>--clean</code></td><td>Flag</td><td><code>False</code></td><td>Clean up the temporary folder after testing.</td></tr></tbody></table>

**Examples:**

**Basic Testing (local model):**

```bash
$ ersilia test eosxxxx --from_dir /path/to/model
```

**Testing with Different Sources:**

```bash
$ ersilia test eosxxxx --from_github/--from_dockerhub/--from_s3
```

**Testing with Specific Levels:**

```bash
$ ersilia test eosxxxx --surface/--shallow/--deep
```
