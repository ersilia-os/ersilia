import json
import tempfile

from ... import ModelBase
from ...core.session import Session
from ...io.input import ExampleGenerator
from .. echo import echo

"""Create example command"""
def example(model, file_name, simple, random, n_samples = 5, deterministic = False):
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

    #if model is not None:
    #     model_id = ModelBase(model).model_id
    # else:
    #     session = Session(config_json=None)
    #     model_id = session.current_model_id()
    # if not model_id:
    #     raise RuntimeError("No model found. Please specify a model or serve a model in the current shell.", fg="red",)
    # eg = ExampleGenerator(model_id=model_id)
    # with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".json") as output:
    #     output_file_name = output.name
    #     eg.example(
    #           n_samples,
    #           output_file_name,
    #           try_predefined=not random,
    #           deterministic=deterministic,
    #           simple=True
    #       )
    #     output.seek(0)
    #     example_data = json.load(output)
    # return example_data