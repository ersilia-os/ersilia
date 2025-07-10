import json
import tempfile
import sys
import os

import pandas as pd

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.exceptions_utils.api_exceptions import UnprocessableInputError
from ...utils.terminal import is_quoted_list
from ..echo import echo


def validate_input_output_types(input, output):
    if not (isinstance(input, list)) or (isinstance(input, str) and input.endswith(".csv")): 
        echo(
            "Input format invalid. Please provide a string, list, and or .csv input instead.",
            fg="red",
            bold=True,
        )
        sys.exit(1)
    if output is not None and not any(
        [output.endswith(ext) for ext in (".csv", ".h5", ".json")]
    ):
        echo(
            "This output type is not allowed in Ersilia. A valid output types are .csv, .h5 or .json",
            fg="red",
            bold=True,
        )
        sys.exit(1)
    # if output is None:
    #     echo(
    #         "Please specify a valid output types which are .csv, .h5 or .json",
    #         fg="red",
    #         bold=True,
    #     )
    #     sys.exit(1)


# input is list
def run(model_id, input, output, batch_size=100):
    # Runs the current model on a list of SMILES strings and
    # returns the prediction as a pandas data frame.

    #   Args:
    # 		input - a list of SMILES strings
    # 		batch_size - number of SMILES to process per batch

    # 	Returns:
    # 		A pandas df with the predictions
    validate_input_output_types(input,output)
    session = Session(config_json=None)
    model_id = session.current_model_id()
    service_class = session.current_service_class()
    if model_id is None:
        print(
            "No model seems to be served. Please run 'ersilia serve ...' before.",
            fg="red",
        )
        return
    # Input is a string or list of SMILEs
    if type(input) == str or isinstance(input, list):
        input_df = pd.DataFrame({"input": input})
        input_file = tempfile.NamedTemporaryFile(mode="w+t", encoding="utf-8", suffix=".csv", delete=False)
        input_df.to_csv(input_file.name, index=False)
        input_file.flush()
    if output is None:
        output_path = os.path.join(os.getcwd(), "output_results.csv")
    else: 
        output_path = str(output)
    mdl = ErsiliaModel(
                model_id,
                output_source=output_path,
                service_class=service_class,
                config_json=None,
    )
    if type(input) == str or isinstance(input, list):
        os.remove(input_file)
    result = mdl.run(input=input_file.name, output=output_path, batch_size=batch_size)

    print(
        f"✅ The output successfully generated in {output_path} file!",
    )

    # output_path.seek(0)

    # Return as dataframe
    # return pd.read_csv(output)
    #   return pd.read_csv(output_file.name)
    # Input is a file
    # if isinstance(input, (str, os.PathLike)) and str(input).lower().endswith("csv") and os.path.isfile(obj):
        # check formatting?