import os
import sys

import rich_click as click

from ... import throw_ersilia_exception
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
    @click.option(
        "--n_samples",
        "-n",
        default=None,
        type=click.IntRange(min=1),
        help="Number of examples to generate. Ignored in curated mode.",
    )
    @click.option(
        "--output_file", "-o", default=None, type=click.STRING, help="Output file name."
    )
    @click.option("--file_name", "-f", default=None, type=click.STRING, hidden=True)
    @click.option(
        "--mode",
        "-m",
        type=click.Choice(["random", "curated", "deterministic"], case_sensitive=False),
        default="random",
        show_default=True,
        help="How examples are generated: random, curated, or deterministic.",
    )
    @throw_ersilia_exception()
    def example(model, n_samples, output_file, file_name, mode):
        from ... import ModelBase
        from ...core.session import Session
        from ...io.input import ExampleGenerator

        # support legacy --file_name / -f
        resolved_file = output_file or file_name
        if not resolved_file:
            raise click.UsageError("Missing option '--output_file' / '-o'.")
        from ..run_checks import fail

        if not resolved_file.lower().endswith(".csv"):
            fail(
                "The output file must end in .csv.",
                "Examples are written as CSV, so they can be used with 'ersilia run'.",
            )
        folder = os.path.dirname(resolved_file) or "."
        if not os.path.isdir(folder):
            fail(f"The output folder {folder} does not exist.")
        if mode == "curated" and n_samples is not None:
            echo("--n_samples is ignored in curated mode.", fg="yellow")
        if n_samples is None and mode != "curated":
            n_samples = 5
        if model is not None:
            model_id = ModelBase(model).model_id
        else:
            session = Session(config_json=None)
            model_id = session.current_model_id()
        if not model_id:
            echo(
                "No model was given, and no model is being served in this terminal.",
                fg="red",
            )
            echo("Give one, e.g. 'ersilia example eos42ez -o input.csv'.")
            sys.exit(1)
        if mode == "curated":
            from ...default import PREDEFINED_EXAMPLE_FILES

            model_dir = ModelBase(model_id)._model_path(model_id)
            if not any(
                os.path.exists(os.path.join(model_dir, f))
                for f in PREDEFINED_EXAMPLE_FILES
            ):
                fail(
                    f"Model {model_id} has no curated examples here.",
                    "Fetch the model first, or use '--mode random' instead.",
                )
        eg = ExampleGenerator(model_id=model_id)
        eg.example(
            n_samples,
            resolved_file,
            mode=mode,
        )
        echo(f"Examples written to {resolved_file}.", fg="green")

    return example
