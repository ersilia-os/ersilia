# Ersilia API
Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called ErsiliaAPI.

To use the API, create a file or open a jupyter notebook. 

# Import the Class
<pre> ```python from ersilia.api import ErsiliaModel ``` </pre>
hello test
# Instantiate the model(ex: Retrosynthetic Accessibility Score)
mdl_retro = ErsiliaModel("eos2r5a")

# Fetch model
mdl_retro.fetch(verbose=False)

# Serve model
mdl_retro.serve(verbose=False)

# Check Fetched Status
mdl_retro.is_fetched()

# Check Docker Status
mdl_retro.is_docker()


# Run Model
input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
]
mdl_retro.run(input, batch_size=100) 

# Example Command
example = molecular_weight.example(True, True, 10, False)
 

# Info Command
mdl_retro.info()


# Close Model
mdl_retro.close()


# Delete Model
mdl_retro.delete()


# With statement: Automatic Serving and Closing
with mdl_retro as model:
    model.info()
    model.run(input, batch_size=100)


# Ersilia Hub Class/Catolog Command
from ersilia.api import ErsiliaHub

Hub = ErsiliaHub()
df = Hub.catalog()
#to convert data frame to csv file
df.to_csv(catalog_file.csv)