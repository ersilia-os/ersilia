import json
import os
import sys
import types

import click

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.exceptions_utils.api_exceptions import UnprocessableInputError
from ...utils.terminal import is_quoted_list
from .. import echo
from . import ersilia_cli


def run_cmd():
    """
    Runs a specified model.

    This command allows users to run a specified model with given inputs.

    Returns
    -------
    function
        The run command function to be used by the CLI and for testing in the pytest.

    Examples
    --------
    .. code-block:: console

        Run a model by its ID with input data:
        $ ersilia run -i <input_data> --as-table

        Run a model with batch size:
        $ ersilia run -i <input_data> -b 50
    """

    def validate_input_output_types(input, output):
        if (type(input) == str and not input.endswith(".csv")) or is_quoted_list(
            json.dumps(input)
        ):
            echo(
                "String and list input types are not allowed in Ersilia. Please a csv input instead",
                fg="red",
                bold=True,
            )
            sys.exit(1)
        if output is not None and not any(
            [output.endswith(ext) for ext in (".csv", ".h5", ".json")]
        ):
            echo(
                "This output type is not allowed in Ersilia. A valid output types are .csv, .h5 or .json",
                fg="red",
                bold=True,
            )
            sys.exit(1)
        if output is None:
            echo(
                "Please specify a valid output types which are .csv, .h5 or .json",
                fg="red",
                bold=True,
            )
            sys.exit(1)

    # Example usage: ersilia run -i {INPUT} [-o {OUTPUT} -b {BATCH_SIZE}]
    @ersilia_cli.command(short_help="Run a served model", help="Run a served model")
    @click.option("-i", "--input", "input", required=True, type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    @click.option(
        "-b", "--batch_size", "batch_size", required=False, default=100, type=click.INT
    )
    def run(input, output, batch_size):
        validate_input_output_types(input, output)
        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        output_source = session.current_output_source()
        if model_id is None:
            echo(
                "No model seems to be served. Please run 'ersilia serve ...' before.",
                fg="red",
            )
            return

        mdl = ErsiliaModel(
            model_id,
            output_source=output_source,
            service_class=service_class,
            config_json=None,
        )
        try:
            result = mdl.run(input=input, output=output, batch_size=batch_size)
            iter_values = []
            if isinstance(result, types.GeneratorType):
                for result in mdl.run(
                    input=input, output=output, batch_size=batch_size
                ):
                    if result is not None:
                        iter_values.append(result)
            echo(
                f"✅ The output successfully generated in {output} file!",
                fg="green",
                bold=True,
            )
        except UnprocessableInputError as e:
            echo(f"❌ Error: {e.message}", fg="red")
            echo(f"💡 {e.hints}")
            if output and os.path.exists(output):
                os.remove(output)
            return

    return run
