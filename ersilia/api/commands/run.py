import os
import sys
import tempfile
import pandas as pd
import types

from ... import ErsiliaModel
from ...core.session import Session
from ..echo import echo


def validate_input_output_types(input, output):
    """
    Validates that 'input' is either a Python list or a path to a .csv file,
    and that 'output' (if provided) ends with a valid extension.
    """
    if not (isinstance(input, list) or (isinstance(input, str) and input.lower().endswith(".csv"))):
        echo(
            "Input format invalid. Please provide a list of SMILEs or a .csv path.",
            fg="red",
            bold=True,
        )
        sys.exit(1)
    
    if output is not None and not any(
        output.lower().endswith(ext) for ext in (".csv", ".h5", ".json")
    ):
        echo(
            "This output type is not allowed in Ersilia. A valid output types are .csv, .h5 or .json",
            fg="red",
            bold=True,
        )
        sys.exit(1)


def run(model_id, input, output, batch_size=100):
    """
    Runs the current model on a list of SMILES strings and
    returns the prediction as a pandas data frame.

    Args
    ----
    input: a list or a path to a CSV file containing SMILES strings.
    batch_size: number of SMILES to process per batch

    Returns
    -------
    Str:    The path to the output file if the output successfully generated..
    #     The run command function to be used by the API.
    #     A pandas df with the predictions.

    """
    validate_input_output_types(input, output)
    session = Session(config_json=None)
    model_id = model_id
    service_class = session.current_service_class()
    output_source = output or session.current_output_source() or os.path.join(os.getcwd(), "output_results.csv")
    
    if model_id is None:
        echo(
            "No model seems to be served. Please run 'ersilia serve ...' before.",
            fg="red",
        )
        return
    
    cleanup_input = False
    if isinstance(input, list):
        # Write list to a temporary CSV
        cleanup_input = True
        df_input = pd.DataFrame({"input": input})
        temp_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
        pd.DataFrame({"input": input}).to_csv(temp_file.name, index=False)
        input_path = temp_file.name
    else:
        # already a CSV file
        input_path = input

    output_path = output or os.path.join(os.getcwd(), "output_results.csv")
    
    mdl = ErsiliaModel(
        model_id,
        output_source=output_path,
        service_class=service_class,
        config_json=None,
    )
    mdl.run(input=input_path, output=output_source, batch_size=batch_size)

    if output_source.lower().endswith(".csv") and os.path.exists(output_source):
        df = pd.read_csv(output_source)
    else:
       echo("Output generated but not in CSV format.", fg="yellow")

    if cleanup_input:
        try:
            os.remove(input_path)
        except OSError:
            pass

    echo(f":check_mark_button: The output successfully generated in {output_path} file!", fg="green", bold=True)

    return df