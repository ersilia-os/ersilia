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
    @ersilia_cli.group(
            short_help="Generate sample of Ersilia models or model inputs",
            help="""This command can sample both ersilia models, or inputs for a given or currently running model.\n
For the model input, the number of examples can be specified, as well as a file name.\n
Simple inputs only contain the essential information, while complete inputs contain key and other fields, potentially.\n
For ersilia models, only model identifiers are returned for a given sample size.
""",
    )
    @click.argument("model", required=False, default=None, type=click.STRING)
    @click.option("--n_samples", "-n", default=5, type=click.INT)
    @click.option("--file_name", "-f", default=None, type=click.STRING)
    @click.option("--simple/--complete", "-s/-c", default=True)
    @click.option("--predefined/--random", "-p/-r", default=True)
    def example(model, n_samples, file_name, simple, predefined):
        if model is not None:
            model_id = ModelBase(model).model_id
        else:
            session = Session(config_json=None)
            model_id = session.current_model_id()
        if not model_id:
            echo(
                "No model found. Please specify a model or serve a model in the current shell.", fg="red"
            )
            return
        eg = ExampleGenerator(model_id=model_id)
        if file_name is None:
            echo(
                json.dumps(
                    eg.example(n_samples, file_name, simple, try_predefined=predefined),
                    indent=4,
                )
            )
        else:
            eg.example(n_samples, file_name, simple, try_predefined=predefined)


    @example.command()
    @click.option("--n_samples", "-n", default=5, type=click.INT)
    @click.option("--file_name", "-f", default=None, type=click.STRING)
    def models(n_samples, file_name):
        sampler = ModelSampler(config_json=None)
        sampler.sample(n_samples=n_samples, file_name=file_name)
        if file_name is None:
            echo(json.dumps(sampler.sample(n_samples), indent=4))
        else:
            sampler.sample(n_samples, file_name)