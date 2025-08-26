import os
import tempfile
import types

import pandas as pd

from ... import ErsiliaModel
from ...core.session import Session
from ..echo import echo


def run(model_id, input_list, batch_size=100):
    """
    Runs the current model on input strings, optionall saving the results to a file and always returning
    a DataFarme.

    Parameters
    ----
    model_id : str
        Identifier of the Ersilia model to invoke.
    input: list of str or str
        A list of input strings or a filepath to a CSV file (.csv) containing inputs
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
    if isinstance(input_list, list):
        cleanup_input = True
        temp_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
        pd.DataFrame({"input": input_list}).to_csv(temp_file.name, index=False)
        input_path = temp_file.name
    else:
        raise Exception("Input must be a list")

    # Prepare output filepath
    cleanup_output = True
    temp_file_output = tempfile.NamedTemporaryFile(
        mode="w+", suffix=".csv", delete=False
    )
    output_path = temp_file_output.name

    # Run the model
    mdl = ErsiliaModel(
        model_id,
        output_source=output_source,
        service_class=service_class,
        config_json=None,
    )
    result = mdl.run(input=input_path, output=output_path, batch_size=batch_size)
    iter_values = []
    if isinstance(result, types.GeneratorType):
        for result in mdl.run(
            input=input_list, output=output_path, batch_size=batch_size
        ):
            if result is not None:
                iter_values.append(result)

    # Notify the user on persistent output
    if not cleanup_output:
        echo(
            f"✅ The output was successfully generated at {output_path}!",
            fg="green",
            bold=True,
        )
    # Load results into DataFrame
    df = pd.read_csv(output_path)

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
