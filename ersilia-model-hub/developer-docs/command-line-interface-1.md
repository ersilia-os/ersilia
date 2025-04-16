---
description: Command Line Interface to interact with the Ersilia Model Hub
---

# Command line interface

The table details all the commands availabe in Ersilia's CLI.

### 1. Fetch Command

This command fetches a specified model from various sources (e.g., GitHub, DockerHub, local directory, AWS S3).

<table><thead><tr><th width="256">Option / Argument</th><th width="102">Type</th><th width="106">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><em>Required</em></td><td>The model ID to be fetched.</td></tr><tr><td><code>--overwrite/--reuse</code></td><td>Boolean</td><td><code>True</code></td><td>Overwrite the existing environment or reuse an already available one for this model.</td></tr><tr><td><code>--from_dir</code></td><td>String</td><td><code>None</code></td><td>Specifies a local path where the model is stored.</td></tr><tr><td><code>--from_github</code></td><td>Flag</td><td><code>False</code></td><td>Fetch directly from GitHub.</td></tr><tr><td><code>--from_dockerhub</code></td><td>Flag</td><td><code>False</code></td><td>Force fetch from DockerHub.</td></tr><tr><td><code>--version</code></td><td>String</td><td><code>None</code></td><td>Version of the model to fetch, used when fetching from DockerHub.</td></tr><tr><td><code>--from_s3</code></td><td>Flag</td><td><code>False</code></td><td>Force fetch from an AWS S3</td></tr></tbody></table>

**Fetch Example:**

```bash
$ ersilia fetch eosxxxx --from_github/--from_dockerhub/--from_s3
```

***

### 2. Serve Command

This command serves a specified model as an API.

<table><thead><tr><th width="237">Option / Argument</th><th width="96">Type</th><th width="113">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><em>Required</em></td><td>The model ID to be served.</td></tr><tr><td><code>--port</code> / <code>-p</code></td><td>Integer</td><td><code>None</code></td><td>Preferred port to use for serving the model.</td></tr><tr><td><code>--track</code> / <code>-t</code></td><td>Flag</td><td><code>False</code></td><td>Enables tracking of the serve session.</td></tr><tr><td><code>--cache/--no-cache</code></td><td>Flag</td><td><code>True</code></td><td>Toggle caching on or off.</td></tr><tr><td><code>--max-cache-memory-frac</code></td><td>Float</td><td><code>None</code></td><td>Sets the maximum fraction of memory to use for caching.</td></tr></tbody></table>

**Serve Example:**

Serve a model on port 8080 and enable tracking:

```bash
$ ersilia serve eosxxxx --no-cache
```

***

### 3. Run Command

This command runs a specified model with given inputs.

<table><thead><tr><th width="189">Option / Argument</th><th width="97">Type</th><th width="118">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>-i, --input</code></td><td>String</td><td><em>Required</em></td><td>The input data for running the model.</td></tr><tr><td><code>-o, --output</code></td><td>String</td><td><code>None</code></td><td>Destination for model output (if applicable).</td></tr><tr><td><code>-b, --batch_size</code></td><td>Integer</td><td><code>100</code></td><td>Batch size to use when running the model.</td></tr><tr><td><code>--as_table</code> / <code>-t</code></td><td>Flag</td><td><code>False</code></td><td>Displays the output as a table format.</td></tr></tbody></table>

***

**Run Example:**

Run a model using specified input data and a batch size of 50, displaying the output as a table:

```bash
$ ersilia run -i input.csv -o output.csv -b 50 --as_table
```

### 4. Example Command

This command generates input examples to be tested for a specific model. It allows you to specify the number of examples, the output file name, and whether the examples are simple or complete, as well as random or predefined.

<table><thead><tr><th width="242">Option / Argument</th><th width="110">Type</th><th width="137">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><code>None</code></td><td>Optional model ID for which to generate examples.</td></tr><tr><td><code>--n_samples</code>, <code>-n</code></td><td>Integer</td><td><code>5</code></td><td>Number of input examples to generate.</td></tr><tr><td><code>--file_name</code>, <code>-f</code></td><td>String</td><td><code>None</code></td><td>File name where the examples should be saved.</td></tr><tr><td><code>--simple</code> / <code>--complete</code></td><td>Flag</td><td><code>True</code> (simple)</td><td>Chooses between simple (essential information only) and complete (includes key and additional fields) examples.</td></tr><tr><td><code>--random</code> / <code>--predefined</code></td><td>Flag</td><td><code>True</code></td><td>Determines if the examples are randomly generated or predefined.</td></tr></tbody></table>

**Example:**

```bash
$ ersilia example eosxxxx -n 10 --file_name examples.csv --complete --random
```

***

### 5. Catalog Command

This command lists a catalog of models available either locally or in the model hub. You can display additional details, view model cards, or choose the output format (table or JSON).

<table><thead><tr><th width="244">Option / Argument</th><th width="78">Type</th><th width="108">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>--hub/--local</code></td><td>Flag</td><td><code>False</code></td><td>Toggle between showing models available in the model hub or locally.</td></tr><tr><td><code>--file_name</code>, <code>-f</code></td><td>String</td><td><code>None</code></td><td>Name of the catalog file to be used.</td></tr><tr><td><code>--more</code> / <code>--less</code></td><td>Flag</td><td><code>False</code></td><td>Display more detailed information than just the EOS identifier.</td></tr><tr><td><code>--card</code></td><td>Flag</td><td><code>False</code></td><td>When provided, displays the model card for the specified model ID.</td></tr><tr><td><code>model</code> (argument)</td><td>String</td><td><em>Optional</em></td><td>Model ID for which to display the model card when using the <code>--card</code> flag.</td></tr><tr><td><code>--as-json /--as-table</code></td><td>Flag</td><td><code>False</code></td><td>Toggle output format between JSON and table format.</td></tr></tbody></table>

**Example:**

```bash
$ ersilia catalog --card eosxxxx --as-json
```

### 6. Delete Command

This command deletes a specified model from the local storage. You can delete a specific model by providing its ID or remove all locally available models using the `--all` flag.

<table><thead><tr><th width="174">Option / Argument</th><th width="74">Type</th><th width="118">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><em>Optional</em></td><td>The model ID to delete.</td></tr><tr><td><code>--all</code></td><td>Flag</td><td><code>False</code></td><td>When provided, deletes all locally available models.</td></tr></tbody></table>

**Examples:**

*   **Delete a specific model:**

    ```bash
    $ ersilia delete eosxxxx
    ```
*   **Delete all models:**

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

### 8. Test Command

This command tests a model and produces performance metrics, allowing you to specify various options such as the source from which to fetch the model, and the depth of the testing process.

<table><thead><tr><th width="182">Option / Argument</th><th width="80">Type</th><th width="109">Default Value</th><th>Description</th></tr></thead><tbody><tr><td><code>model</code> (argument)</td><td>String</td><td><em>Required</em></td><td>The model ID to test.</td></tr><tr><td><code>--from_dir</code></td><td>String</td><td><code>None</code></td><td>Local directory where the model is stored.</td></tr><tr><td><code>--from_github</code></td><td>Flag</td><td><code>False</code></td><td>Fetch the model directly from GitHub.</td></tr><tr><td><code>--from_dockerhub</code></td><td>Flag</td><td><code>False</code></td><td>Force fetching the model from DockerHub.</td></tr><tr><td><code>--from_s3</code></td><td>Flag</td><td><code>False</code></td><td>Force fetching the model from an AWS S3 bucket.</td></tr><tr><td><code>--version</code></td><td>String</td><td><code>None</code></td><td>Specify the version of the model to fetch, especially useful when fetching from DockerHub.</td></tr><tr><td><code>--shallow</code></td><td>Flag</td><td><code>False</code></td><td>Run shallow tests (e.g., container size, output consistency).</td></tr><tr><td><code>--deep</code></td><td>Flag</td><td><code>False</code></td><td>Run deep tests (e.g., computational performance checks).</td></tr><tr><td><code>--surface</code></td><td>Flag</td><td><code>False</code></td><td>Additional flag for testing, similar to deep checks.</td></tr><tr><td><code>--inspect</code></td><td>Flag</td><td><code>False</code></td><td>Run inspection tests for further model diagnostics.</td></tr><tr><td><code>--report_path</code></td><td>String</td><td><code>None</code></td><td>Specify the file path to output the test report in JSON format.</td></tr><tr><td><code>--clean</code></td><td>Flag</td><td><code>False</code></td><td>Clean up the temporary folder after testing.</td></tr></tbody></table>

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
