# API
Models can be fetched from the Ersilia Model Hub, served, and run as a Python package. The main class is called ErsiliaAPI:

# import main class
from ersilia.api import ErsiliaModel
# instantiate the model(ex: Retrosynthetic Accessibility Score)
# name the model(ex : mdl_retro)
mdl_retro = ErsiliaModel("eos2r5a")
Then, you can perform the same actions as in the CLI. Specify preference for verbosity of status messages (automatically succinct). To fetch and serve:

# fetch model
mdl_retro.fetch(verbose=False)

#serve model
mdl_retro.serve(verbose=False)
To check if a model is fetched, use is_fetched:

mdl_retro.is_fetched()
To check if a docker is running locally, use is_docker:

mdl_retro.is_docker()
To make predictions for Halicin and Ibuprofen:

# Halicin and Ibuprofen
input = [
    "C1=C(SC(=N1)SC2=NN=C(S2)N)[N+](=O)[O-]",
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
]
# predict
mdl_retro.run(input, batch_size=100)
# Pass the input for the model as a .csv file with only one column with header or a list of SMILES compounds
# Flexibility to specify file path for desired output file or automatically create an output file in the current directory
# Specify the batch size for generating model predictions. By default, Ersilia works with batch size of 100 inputs.
To sample inputs for a given model, use the example command. 


example = molecular_weight.example(True, True, 10, False)
#example command returns a data frame
#to convert to a csv
example.to_csv(example_file.csv) 
# Parameters: 
# mdl.example(file_name, simple, random, n_samples, deterministic)
# file_name: str (File name where the examples should be saved.)
# simple: bool (Simple inputs only contain the SMILES, while complete inputs also include InChIKey and the molecule's name.)
# random: bool(If the model source contains an example input file, when the predefined flag is set, then inputs are sampled from that file. Only the number of samples present in the file are returned, especially if --n_samples is greater than that number. By default, Ersilia samples inputs randomly.
# n_samples: int(Specify the number of example inputs to generate for the given model.)
# deterministic: bool (Used to generate examples data deterministically instead of random sampling. This allows when every time you run with example command with this flag you get the same types of examples.)
To get detailed information about a current active session, use the info command. 


mdl_retro.info()
To close the model:


# close model
mdl_retro.close()
To delete the model:


# delete model
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

# use with statement
# this allows for automatic serving and closing of an already fetched model
with mdl_retro as model:
    model.info()
    model.run(input, batch_size=100)
Using the catalog command
This command allows users to access a catalog of models available either locally or in the model hub.


# to use this command, import the ErsiliaHub class. 
from ersilia.api import ErsiliaHub

Hub = ErsiliaHub()
df = Hub.catalog()
#to convert data frame to csv file
df.to_csv(catalog_file.csv)