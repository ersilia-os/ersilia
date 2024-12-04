import os
import click
import json
import tempfile

from ...cli import echo
from . import ersilia_cli
from ersilia.cli.commands.run import run_cmd
from ersilia.core.base import ErsiliaBase
from ...publish.test import ModelTester

from ersilia.utils.exceptions_utils import throw_ersilia_exception

from ersilia.utils.exceptions_utils.test_exceptions import WrongCardIdentifierError

from ersilia.default import INFORMATION_FILE


def test_cmd():
    """Test a model"""

    # Example usage: ersilia test {MODEL}
    @ersilia_cli.command(
        short_help="Test a model",
        help="Test a model and obtain performance metrics",
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "-o", "--output", "output", required=False, default=None, type=click.STRING
    )
    def test(model, output):
        mdl = ModelTester(model)
        model_id = mdl.model_id

        if model_id is None:
            echo(
                "No model seems to be served. Please run 'ersilia serve ...' before.",
                fg="red",
            )
            return

        mt = ModelTester(model_id=model_id)
        # click.echo("Checking model information")
        mt.run(output)  # pass in the output here
