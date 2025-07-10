import json

from ... import ModelBase
from ...core.session import Session
from ...io.input import ExampleGenerator
from ..echo import echo


def example(model, file_name, simple=True, random=True, n_samples=5, deterministic=False):  
    """
    This command can sample inputs for a given model. 
    
    Args
    -------
    model: The model ID to be served. Can either be the eos identifier or the slug identifier.
    file_name: File name where the examples should be saved.
    simple: Simple inputs only contain the SMILES, while complete inputs also include InChIKey and the molecule's name.
    random: If the model source contains an example input file, when the predefined flag is set, then inputs are sampled from that file. Only the number of samples present in the file are returned, especially if --n_samples is greater than that number. By default, Ersilia samples inputs randomly.
    n_samples: Specify the number of example inputs to generate for the given model.
    deterministic: Used to generate examples data deterministically instead of random sampling. This allows when every time you run with example command with this flag you get the same types of examples. 

    Returns
    -------
    Function: The fetch command function to be used by the CLI and for testing in the pytest.
    Str: Confirmation message on success or warning message on failure.
    
    Raises
    -------
    RuntimeError: If both BentoML and FastAPI are used togehter. 

    """
    if model is not None:
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
    eg = ExampleGenerator(model_id=model_id)
    if file_name is None:
        echo(
            json.dumps(
                eg.example(
                    n_samples,
                    file_name,
                    simple,
                    try_predefined=not random,
                    deterministic=deterministic,
                ),
                indent=4,
            )
        )
    else:
        eg.example(
            n_samples,
            file_name,
            simple,
            try_predefined=not random,
            deterministic=deterministic,
        )
    return example

    