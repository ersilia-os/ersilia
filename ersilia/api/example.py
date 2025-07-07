import json
import tempfile

from ... import ModelBase
from ..core.session import Session
from ..io.input import ExampleGenerator
from .. import echo
from . import ersilia_cli

"""Create example command"""
def example(model, n_samples=5, random=True, deterministic=True):
    if model is not None:
        model_id = ModelBase(model).model_id
    else:
        session = Session(config_json=None)
        model_id = session.current_model_id()
    if not model_id:
        raise RuntimeError("No model found. Please specify a model or serve a model in the current shell.", fg="red",)
    eg = ExampleGenerator(model_id=model_id)
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".json") as output:
        output_file_name = output.name
        eg.example(
              n_samples,
              output_file_name,
              try_predefined=not random,
              deterministic=deterministic,
          )
        output.seek(0)
        example_data = json.load(output)
    return example_data