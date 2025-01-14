import json
import types

import click

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.terminal import print_result_table
from .. import echo
from . import ersilia_cli


def truncate_output(output, max_length=10):
    """
    Truncates long outputs (like arrays or lists) for better readability.
    
    Parameters
    ----------
    output : Any
        The result to be truncated if it is a list or array-like object.
    max_length : int, optional
        The maximum number of elements to display, by default 10.
    
    Returns
    -------
    str
        The truncated output as a formatted string.
    """
    if isinstance(output, (list, tuple)) and len(output) > max_length:
        truncated = output[:max_length]
        return json.dumps(truncated, indent=4) + f"\n... (and {len(output) - max_length} more elements)"
    elif isinstance(output, dict):
        truncated_dict = {k: truncate_output(v, max_length) for k, v in output.items()}
        return json.dumps(truncated_dict, indent=4)
    else:
        return json.dumps(output, indent=4)


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

    # Example usage: ersilia run -i {INPUT} [-o {OUTPUT} -b {BATCH_SIZE}]
    @ersilia_cli.command(short_help="Run a served model", help="Run a served model")
    @click.option("-i", "--input", "input", required=True, type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    @click.option(
        "-b", "--batch_size", "batch_size", required=False, default=100, type=click.INT
    )
    @click.option("--as-table/-t", is_flag=True, default=False)
    def run(input, output, batch_size, as_table):
        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        track_runs = session.tracking_status()

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
            track_runs=track_runs,
        )
        result = mdl.run(
            input=input,
            output=output,
            batch_size=batch_size,
            track_run=track_runs,
        )
        if isinstance(result, types.GeneratorType):
            for result in mdl.run(input=input, output=output, batch_size=batch_size):
                if result is not None:
                    if as_table:
                        print_result_table(result)
                    else:
                        formatted = truncate_output(result)
                        echo(formatted)
                else:
                    echo("Something went wrong", fg="red")
        else:
            if as_table:
                print_result_table(result)
            else:
                try:
                    truncated = truncate_output(result)
                    echo(truncated)
                except Exception as e:
                    echo(f"Error formatting result: {e}", fg="red")

    return run
