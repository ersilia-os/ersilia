# Ersilia API

Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called ErsiliaAPIModel.

To use the API, create a file or open a jupyter notebook. 

## Import the Class
<pre> from ersilia.api.create_api import ErsiliaAPIModel </pre>

Then, you can perform the same actions as in the CLI. Specify preference for the verbosity of status messages (automatically succinct). 

## Instantiate the model(ex: Retrosynthetic Accessibility Score)
<pre>
mdl_retro = ErsiliaAPIModel('eos2r5a')
</pre>

## Fetch model
<pre>
mdl_retro.fetch(verbose=False)
</pre>

## Serve model
<pre>
mdl_retro.serve(verbose=False)
</pre>

## Check Fetched Status
<pre>
mdl_retro.is_fetched()
</pre>

## Check Docker Status
To check if a docker is running locally, use is_docker command:
<pre> mdl_retro.is_docker() </pre>

## Run Model
<pre>
input = [
    'C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]',
    'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O'
]
df = mdl_retro.run(input, None, batch_size=100)
print(df)
</pre>

Alternatively, you can input a CSV file instead of a list. 
<pre>
mdl_retro.run('molecules.csv', 'output.csv', 100)
</pre>

## Example Command
<pre>
example = molecular_weight.example(True, True, 10, False)
</pre>

## Info Command
<pre>
mdl_retro.info()
</pre>

## Close Model
<pre>
mdl_retro.close()
</pre>

## Delete Model
<pre>
mdl_retro.delete()
</pre>

## With statement: Automatic Serving and Closing
A more concise way to run models would be to use the with clause:
<pre>
with mdl_retro as model:
    model.info()
    model.run(input, 'output.csv', batch_size=100)
</pre>

## Ersilia Hub Class/Catalog Command
This command allows users to access a catalog of models available either locally or in the model hub. 
In order to use the catalog command, import the class Ersilia Catalog. 
<pre>
from ersilia.api import ErsiliaCatalog

Hub = ErsiliaCatalog()
df = Hub.catalog()  # Default: shows local catalog
df = Hub.catalog(hub=True)  # Shows hub catalog
df = Hub.catalog(hub=False)  # Shows local catalog (default)
</pre>