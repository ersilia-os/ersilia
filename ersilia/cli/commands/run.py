import sys

import click

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

    # Example usage: ersilia run -i {INPUT} [-o {OUTPUT} -b {BATCH_SIZE}]
    @ersilia_cli.command(
        short_help="Run predictions on the served model",
        help="Run predictions using the currently served model. Input must be a single-column CSV file. Output can be saved as .csv or .h5. A model must be served before running this command.",
    )
    @click.option(
        "-i",
        "--input",
        "input",
        required=True,
        type=click.STRING,
        help="Path to a single-column CSV file containing the input data.",
    )
    @click.option(
        "-o",
        "--output",
        "output",
        required=True,
        default=None,
        type=click.STRING,
        help="Path to the output file. Accepted formats: .csv, .h5.",
    )
    @click.option(
        "-b",
        "--batch_size",
        "batch_size",
        required=False,
        default=100,
        type=click.IntRange(min=1),
        help="Number of inputs processed per batch.",
    )
    def run(input, output, batch_size):
        import os
        import re

        from ... import ErsiliaModel
        from ...core.session import Session
        from ..run_checks import check_run_arguments

        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        output_source = session.current_output_source()

        if model_id is None:
            from ..messages import no_model_served

            no_model_served()
            sys.exit(1)

        output_basename = os.path.basename(output)
        output_model_ids = re.findall(r"eos[0-9][a-z0-9]{3}", output_basename)
        if output_model_ids and output_model_ids[0] != model_id:
            echo(
                f"The output file name mentions {output_model_ids[0]}, but the served model is {model_id}.",
                fg="red",
            )
            echo("Use an output file name that matches the served model.")
            sys.exit(1)

        run_input, tmp_dir = check_run_arguments(input, output)
        mdl = ErsiliaModel(
            model_id,
            output_source=output_source,
            service_class=service_class,
            config_json=None,
        )
        import time

        started = time.time()
        try:
            mdl.run(input=run_input, output=output, batch_size=batch_size)
        finally:
            if tmp_dir:
                import shutil

                shutil.rmtree(tmp_dir, ignore_errors=True)
        # Success means an output file written by this run, not an old one.
        if not os.path.isfile(output) or os.path.getmtime(output) < started - 1:
            echo("No output was written.", fg="red")
            sys.exit(1)
        echo(f"Output written to {output}.", fg="green")

    return run
