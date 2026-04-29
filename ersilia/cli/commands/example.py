import rich_click as click

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
        help="Generate input examples for a model. A model must be served before running this command, or a model identifier must be provided. The number of examples and output file can be specified.",
    )
    @click.argument("model", required=False, default=None, type=click.STRING)
    @click.option("--n_samples", "-n", default=None, type=click.INT, help="Number of examples to generate. Ignored in curated mode.")
    @click.option("--output_file", "-o", default=None, type=click.STRING, help="Output file name.")
    @click.option("--file_name", "-f", default=None, type=click.STRING, hidden=True)
    @click.option(
        "--mode",
        "-m",
        type=click.Choice(
            ["random", "curated", "deterministic"], case_sensitive=False
        ),
        default="random",
        show_default=True,
        help="How examples are generated: random, curated, or deterministic.",
    )
    def example(model, n_samples, output_file, file_name, mode):
        # support legacy --file_name / -f
        resolved_file = output_file or file_name
        if not resolved_file:
            raise click.UsageError("Missing option '--output_file' / '-o'.")
        if mode == "curated" and n_samples is not None:
            echo("Warning: --n_samples is ignored in curated mode.", fg="yellow")
        if n_samples is None and mode != "curated":
            n_samples = 5
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
            resolved_file,
            mode=mode,
        )

    return example
