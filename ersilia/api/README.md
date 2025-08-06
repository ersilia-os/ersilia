# Ersilia API

Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called ErsiliaAPIModel.

To use the API, create a file or open a jupyter notebook. 

## Import the Class
<pre> from ersilia.api import ErsiliaAPIModel ``` </pre>

## Instantiate the model(ex: Retrosynthetic Accessibility Score)
<pre>
mdl_retro = ErsiliaAPIModel("eos2r5a")
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
<pre>
mdl_retro.is_docker()
</pre>

## Run Model
<pre>
input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
]
mdl_retro.run(input, batch_size=100)
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
<pre>
with mdl_retro as model:
    model.info()
    model.run(input, batch_size=100)
</pre>

## Ersilia Hub Class/Catolog Command
<pre>
from ersilia.api import ErsiliaCatalog

Hub = ErsiliaCatalog()
df = Hub.catalog()
</pre>
### to convert data frame to csv file
<pre>
df.to_csv(catalog_file.csv)
</pre>