# imports
import os
import csv
import json
import sys

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

# current file directory
root = os.path.dirname(os.path.abspath(__file__))

# checkpoints directory
checkpoints_dir = os.path.abspath(os.path.join(root, "..", "..", "checkpoints"))

# read checkpoints (here, simply an integer number: 42)
with open(os.path.join(checkpoints_dir, "data.json"), "r") as f:
    ckpt = json.load(f)

# model to be run (repetition of the same JSON output)
def my_model(smiles_list, ckpt):
    return [ckpt for _ in smiles_list]


# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model
outputs = my_model(smiles_list, ckpt)

# write output in a .json file
with open(output_file, "w") as f:
    json.dump(outputs, f, indent=4)
