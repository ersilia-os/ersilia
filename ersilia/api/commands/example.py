import os
import tempfile

import pandas as pd

from ...core.session import Session
from ...io.input import ExampleGenerator
from ..echo import echo


def example(n_samples, mode):
    """
    This command can sample inputs for a given model and save them as a CSV file.

    Args
    -------
    n_samples: Specify the number of example inputs to generate for the given model.
    mode: random, deterministic or predefined
    Returns
    -------
    Function: The exmaple command function to be used by the API.
    Str: Error message if no model was served in the current session.


    """
    session = Session(config_json=None)
    model_id = session.current_model_id()

    if not model_id:
        echo(
            "No model found. Please serve a model in the current shell.",
            fg="red",
        )
        return

    output_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
    eg = ExampleGenerator(model_id=model_id)
    eg.example(n_samples, file_name=output_file.name, mode=mode)
    df = pd.read_csv(output_file.name)
    try:
        output_file.close()
        os.remove(output_file.name)
    except OSError:
        pass
    return df["input"].tolist()
