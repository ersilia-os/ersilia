import click
import json

from . import ersilia_cli
from .. import echo
from ...io.input import ExampleGenerator
from ...db.hubdata.samplers import ModelSampler
from ...core.session import Session
from ... import ModelBase


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
    def example():
        pass

    
    @example.command()
    @click.argument("model", required=False, default=None, type=click.STRING)
    @click.option("--n_samples", "-n", default=5, type=click.INT)
    @click.option("--file_name", "-f", default=None, type=click.STRING)
    @click.option("--simple/--complete", "-s/-c", default=True)
    @click.option("--predefined/--random", "-p/-r", default=True)
    def inputs(model, n_samples, file_name, simple, predefined):
        if model is not None:
            model_id = ModelBase(model).model_id
        else:
            session = Session(config_json=None)
            model_id = session.current_model_id()
            if model_id is None:
                click.echo(
                    click.style("Error: No model id given and no model found running in this shell.", fg="red"),
                    err=True,
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