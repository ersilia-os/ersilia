import os
import sys
import tempfile
import types

import pandas as pd

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.exceptions_utils.api_exceptions import UnprocessableInputError
from ..echo import echo


def validate_input_output_types(input, output):
    """
    Validates that 'input' is either a Python list or a path to a .csv file,
    and that 'output' (if provided) ends with a valid extension.
    """
    if not (isinstance(input, list) or (isinstance(input, str) and input.lower().endswith(".csv"))):
        echo(
            "Input format invalid. Please provide a list of SMILEs or a .csv path.",
            fg="red", bold=True,
        )
        sys.exit(1)
    if (output is not None and not any(output.endswith(ext) for ext in (".csv", ".h5", ".json"))):
        echo(
            "Invalid output type. Valid types are: .csv, .h5, or .json",
            fg="red", bold=True,
        )
        sys.exit(1)

def load_output_to_df(output):
    """
    Loads model output from file to pandas DataFrame based on extension.
    """
    try:
        if output.endswith(".csv"):
            return pd.read_csv(output)
        elif output.endswith(".json"):
            return pd.read_json(output)
        elif output.endswith(".h5"):
            return pd.read_hdf(output)
        else:
            raise ValueError(f"Unsupported output format: {output}")
    except Exception as e:
        echo(f"❌ Failed to load output file {output}: {e}", fg="red")
        raise

def run(model_id, input, output, batch_size=100):
    """
    Runs the current model on a list of SMILES strings and
    returns the prediction as a pandas data frame.

    Args
    ----
    input: a list or a path to a CSV file containing SMILES strings.
    output: path to save the output (must end with .csv, .json, or .h5)
    batch_size: number of SMILES to process per batch

    Returns
    -------
    Str:    The path to the output file if the output successfully generated..
    #     The run command function to be used by the API.
    #     A pandas df with the predictions.

    """
    validate_input_output_types(input, output)
    session = Session(config_json=None)
    service_class = session.current_service_class()
    # output_source = output or session.current_output_source() or os.path.join(os.getcwd(), "output_results.csv")

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
        temp_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
        pd.DataFrame({"input": input}).to_csv(temp_file.name, index=False)
        input_path = temp_file.name
    else:
        # already a CSV file
        input_path = input
    # output_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
    # output_path = output or os.path.join(os.getcwd(), "output_results.csv")

    mdl = ErsiliaModel(
        model_id,
        output_source=output,
        service_class=service_class,
        config_json=None,
    )
    try:
        result = mdl.run(input=input_path, output=output, batch_size=batch_size)
        iter_values = []
        if isinstance(result, types.GeneratorType):
            for result in mdl.run(
                input=input, output=output, batch_size=batch_size
            ):
                if result is not None:
                    iter_values.append(result)
        echo(f"✅ The output was successfully generated at {output}!", fg="green", bold=True)
        # if not os.path.exists(output) or os.path.getsize(output) == 0:
        #     echo(f"❌ Output file {output} is empty or missing.", fg="red")
        #     raise ValueError(f"Output file {output} is empty or missing.")
        df = load_output_to_df(output)
    except UnprocessableInputError as e:
        echo(f"❌ Error: {e.message}", fg="red")
        echo(f"💡 {e.hints}")
        if output and os.path.exists(output):
            os.remove(output)
        raise e

    finally:
        if cleanup_input:
            try:
                os.remove(input_path)
            except OSError:
                pass

    return df
