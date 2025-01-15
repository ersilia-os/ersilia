import json
import types

import click

from ... import ErsiliaModel
from ...core.session import Session
from ...utils.terminal import print_result_table
from .. import echo
from . import ersilia_cli

def truncate_output(output, max_items=10):
    """
    Truncates long outputs for clarity.

    Parameters
    ----------
    output : Any
        The output to process and truncate if necessary.
    max_items : int, optional
        The maximum number of items to display for arrays/lists.

    Returns
    -------
    str
        The truncated output as a formatted string.
    """
    if isinstance(output, list):
        # Truncate lists to the first `max_items` items
        if len(output) > max_items:
            return f"{output[:max_items]} ... (and {len(output) - max_items} more items)"
        return str(output)
    elif isinstance(output, dict):
        # Convert dict to JSON and truncate if necessary
        formatted = json.dumps(output, indent=4)
        lines = formatted.splitlines()
        if len(lines) > max_items:
            return "\n".join(lines[:max_items]) + f"\n... (and {len(lines) - max_items} more lines)"
        return formatted
    elif isinstance(output, str):
        # Truncate strings with a character limit
        if len(output) > 500:
            return f"{output[:500]}... (truncated)"
        return output
    else:
        return str(output)


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

    # Process the result
    if isinstance(result, types.GeneratorType):
        for result in mdl.run(input=input, output=output, batch_size=batch_size):
            if result is not None:
                formatted = truncate_output(result)
                if as_table:
                    print_result_table(formatted)
                else:
                    echo(formatted)
            else:
                echo("Something went wrong", fg="red")
    else:
        # Truncate long outputs
        formatted = truncate_output(result)
        if as_table:
            print_result_table(formatted)
        else:
            try:
                echo(formatted)
            except:
                print_result_table(formatted)

    return run
