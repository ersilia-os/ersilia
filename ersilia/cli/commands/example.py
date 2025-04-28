import json

import click

from ... import ModelBase
from ...core.session import Session
from ...io.input import ExampleGenerator
from .. import echo
from . import ersilia_cli


def example_cmd():
    """Create example command"""

    # Example usage: ersilia example {MODEL} -n 10 [--file_name {FILE_NAME} --simple/--complete]
    @ersilia_cli.command(
        short_help="Generate input examples for the model of interest",
        help="This command generates input examples to be tested for the model of interest. The number of examples can be specified, as well as a file name. Simple inputs only contain the essential information, while complete inputs contain key and other fields, potentially.",
    )
    @click.argument("model", required=False, default=None, type=click.STRING)
    @click.option("--n_samples", "-n", default=5, type=click.INT)
    @click.option("--file_name", "-f", default=None, type=click.STRING)
    @click.option("--simple/--complete", "-s/-c", default=True)
    @click.option("--random/--predefined", "-r/-p", default=True)
    @click.option("--deterministic", "-d", default=False, is_flag=True)
    def example(model, n_samples, file_name, simple, random, deterministic):
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
