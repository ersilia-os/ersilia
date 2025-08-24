import os
import tempfile

import pandas as pd

from ... import ModelBase
from ...core.session import Session
from ...io.input import ExampleGenerator
from ..echo import echo


def example(model, simple=True, random=True, n_samples=5, deterministic=False):
    """
    This command can sample inputs for a given model and save them as a CSV file.

    Args
    -------
    model: The model ID to be served. Can either be the eos identifier or slug.
    file_name: Path where the CSV examples should be saved. Must end with .csv.
    simple: If True, only input strings are returned. If False, outputs include key and input.
    random: If the model source contains an example input file, when the predefined flag is set, then inputs are sampled from that file. Only the number of samples present in the file are returned, especially if --n_samples is greater than that number. By default, Ersilia samples inputs randomly.
    n_samples: Specify the number of example inputs to generate for the given model.
    deterministic: Used to generate examples data deterministically instead of random sampling. This allows when every time you run with example command with this flag you get the same types of examples.

    Returns
    -------
    Function: The exmaple command function to be used by the API.
    Str: Error message if no model was served in the current session.


    """
    if model:
        model_id = ModelBase(model).model_id
    else:
        session = Session(config_json=None)
        model_id = session.current_model_id()

    if not model_id:
        echo(
            "No model found. Please specify a model or serve a model in the current shell.",
            fg="red",
        )
        return

    output_file = tempfile.NamedTemporaryFile(mode="w+", suffix=".csv", delete=False)
    eg = ExampleGenerator(model_id=model_id)
    eg.example(
        n_samples,
        file_name=output_file.name,
        simple=simple,
        try_predefined=not random,
        deterministic=deterministic,
    )
    df = pd.read_csv(output_file.name)
    try:
        output_file.close()
        os.remove(output_file.name)
    except OSError:
        pass
    return df
