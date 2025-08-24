import os
import sys
import tempfile
import types

import pandas as pd

from ... import ErsiliaModel
from ...core.session import Session
from ..echo import echo


def validate_input_output_types(input, output):
    """
    Validate that the provided input and output parameters are of supported types.

    Parameters
    ----------
    input : list of str or str
        A list of input strings or a filepath pointing to a CSV file (.csv) containing inputs with 'input' header
    output : str or None
        Optional filepath for saving the model output. If provided, must end with one of: .csv, .json, .h5

    Raises
    ------
        If 'input' is neither a list nor a .csv filepath, or if 'output' (when not None) does not end with a supported
        extension.
    """
    if not (
        isinstance(input, list)
        or (isinstance(input, str) and input.lower().endswith(".csv"))
    ):
        echo(
            "Input format invalid. Please provide a list of inputs or a .csv path.",
            fg="red",
            bold=True,
        )
        sys.exit(1)
    if output is not None and not any(
        output.endswith(ext) for ext in (".csv", ".h5", ".json")
    ):
        echo(
            "Invalid output type. Valid types are: .csv, .h5, or .json",
            fg="red",
            bold=True,
        )
        sys.exit(1)


def run(model_id, input, output=None, batch_size=100):
    """
    Runs the current model on input strings, optionall saving the results to a file and always returning
    a DataFarme.

    Parameters
    ----
    model_id : str
        Identifier of the Ersilia model to invoke.
    input: list of str or str
        A list of input strings or a filepath to a CSV file (.csv) containing inputs
    output: str, optional
        Filepath to save the model predictions, Supported extensions: .csv, .json, .h5
        If None, a temporary CSV file will be generated and cleaned up after loading.
    batch_size: int, default=100
        Number of input rows to process per batch

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the model's prediction results.

    Raises
    ------
    SystemExit
        If input or output validation fails.

    """
    validate_input_output_types(input, output)
    session = Session(config_json=None)
    service_class = session.current_service_class()
    output_source = session.current_output_source()

    if model_id is None:
        echo(
            "No model seems to be served. Please run 'ersilia serve ...' before.",
            fg="red",
        )
        return

    cleanup_input = False
    cleanup_output = False
    # Prepare input filepath
    if isinstance(input, list):
        cleanup_input = True
        temp_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
        pd.DataFrame({"input": input}).to_csv(temp_file.name, index=False)
        input_path = temp_file.name
    else:
        input_path = input

    # Prepare output filepath
    if output == None:
        cleanup_output = True
        temp_file_output = tempfile.NamedTemporaryFile(
            mode="w+", suffix=".csv", delete=False
        )
        output_path = temp_file_output.name
        output = output_path

    # Run the model
    mdl = ErsiliaModel(
        model_id,
        output_source=output_source,
        service_class=service_class,
        config_json=None,
    )
    result = mdl.run(input=input_path, output=output, batch_size=batch_size)
    iter_values = []
    if isinstance(result, types.GeneratorType):
        for result in mdl.run(input=input, output=output, batch_size=batch_size):
            if result is not None:
                iter_values.append(result)

    # Notify the user on persistent output
    if not cleanup_output:
        echo(
            f"✅ The output was successfully generated at {output}!",
            fg="green",
            bold=True,
        )
    # Load results into DataFrame
    df = pd.read_csv(output)

    # Clean up temporary files
    if cleanup_input:
        try:
            os.remove(input_path)
        except OSError:
            pass
    if cleanup_output:
        try:
            os.remove(output_path)
        except OSError:
            pass

    return df
