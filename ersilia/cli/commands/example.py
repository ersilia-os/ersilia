import click

from ... import ModelBase
from ...core.session import Session
from ...io.input import ExampleGenerator
from .. import echo
from . import ersilia_cli


def example_cmd():
    """Create example command"""

    # Example usage: ersilia example {MODEL} -n 10 [--file_name {FILE_NAME}]
    @ersilia_cli.command(
        short_help="Generate input examples for the model of interest",
        help="This command generates input examples to be tested for the model of interest. The number of examples can be specified, as well as a file name.",
    )
    @click.argument("model", required=False, default=None, type=click.STRING)
    @click.option("--n_samples", "-n", default=5, type=click.INT)
    @click.option("--file_name", "-f", required=True, type=click.STRING)
    @click.option(
        "--mode",
        "-m",
        type=click.Choice(
            ["random", "predefined", "deterministic"], case_sensitive=False
        ),
        default="random",
        show_default=True,
        help="Choose how examples are generated (random, predefined or deterministic).",
    )
    def example(model, n_samples, file_name, mode):
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
        eg.example(
            n_samples,
            file_name,
            mode=mode,
        )

    return example
