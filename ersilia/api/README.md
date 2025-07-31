# API
Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called ErsiliaAPI.

To use the API, create a file or open a jupyter notebook. 

# Import the class
from ersilia.api import ErsiliaModel
# Instantiate the model(ex: Retrosynthetic Accessibility Score)
mdl_retro = ErsiliaModel("eos2r5a")
Then, you can perform the same actions as in the CLI. Specify preference for verbosity of status messages (automatically succinct). To fetch and serve:

# Fetch model
mdl_retro.fetch(verbose=False)

# Serve model
mdl_retro.serve(verbose=False)

# Check Fetched Status

mdl_retro.is_fetched()

# Check Docker Status

mdl_retro.is_docker()

To make predictions for Halicin and Ibuprofen:

input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
]
# predict
mdl_retro.run(input, batch_size=100)
-> Pass the input for the model as a .csv file with only one column with header or a list of SMILES compounds
-> Flexibility to specify file path for desired output file or automatically create an output file in the current directory
-> Specify the batch size for generating model predictions. By default, Ersilia works with batch size of 100 inputs.

To sample inputs for a given model, use the example command. 

# Example Command

example = molecular_weight.example(True, True, 10, False)
#example command returns a data frame
#to convert to a csv
example.to_csv(example_file.csv) 

#Info Command

mdl_retro.info()
To close the model:


# Close Model
mdl_retro.close()
To delete the model:


# Delete Model
mdl_retro.delete()
Using the with statement
A more concise way to run prediction would be to use the with clause:


mdl.fetch()
input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    ]

# As opposed to the following:
mdl.serve()
mdl.info()
mdl.run(input, batch_size=100)
model.close()

# With statement
this allows for automatic serving and closing of an already fetched model
with mdl_retro as model:
    model.info()
    model.run(input, batch_size=100)
Using the catalog command
This command allows users to access a catalog of models available either locally or in the model hub.

# Ersilia Hub Class
to use this command, import the ErsiliaHub class. 
from ersilia.api import ErsiliaHub

Hub = ErsiliaHub()
df = Hub.catalog()
#to convert data frame to csv file
df.to_csv(catalog_file.csv)