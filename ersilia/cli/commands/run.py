import json
import types

import click

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.terminal import print_result_table
from .. import echo
from . import ersilia_cli


def summarize_output(output, max_items=20):
    """
    Summarize the output to display only the first few elements.
    If the output is a list or dictionary, truncate it for better readability.

    Parameters
    ----------
    output : dict, list, or any object
        The output to summarize.
    max_items : int, optional
        The maximum number of items to display, by default 10.

    Returns
    -------
    str
        The summarized output as a string.
    """
    if isinstance(output, list):
        if len(output) > max_items:
            return json.dumps(output[:max_items] + ["..."], indent=4)
    elif isinstance(output, dict):
        if len(output) > max_items:
            truncated = dict(list(output.items())[:max_items])
            truncated["..."] = f"+{len(output) - max_items} more items"
            return json.dumps(truncated, indent=4)
    return json.dumps(output, indent=4)


def run_cmd():
    """
    Runs a specified model.
    """

    @ersilia_cli.command(short_help="Run a served model", help="Run a served model")
    @click.option("-i", "--input", "input", required=True, type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    @click.option(
        "-b", "--batch_size", "batch_size", required=False, default=100, type=click.INT
    )
    @click.option("--as-table/-t", is_flag=True, default=False)
    @click.option("--verbose", is_flag=True, help="Show the full output.")
    def run(input, output, batch_size, as_table, verbose):
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

        def process_result(result):
            if verbose:
                return json.dumps(result, indent=4)
            else:
                return summarize_output(result)

        if isinstance(result, types.GeneratorType):
            for result in mdl.run(input=input, output=output, batch_size=batch_size):
                if result is not None:
                    formatted = process_result(result)
                    if as_table:
                        print_result_table(formatted)
                    else:
                        echo(formatted)
                else:
                    echo("Something went wrong", fg="red")
        else:
            formatted = process_result(result)
            if as_table:
                print_result_table(formatted)
            else:
                try:
                    echo(formatted)
                except Exception:
                    print_result_table(result)

    return run
