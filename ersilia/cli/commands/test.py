import click
from ...cli import echo
from . import ersilia_cli
from ...publish.test import ModelTester


def test_cmd():
    """
    Test a model and obtain performance metrics.

    This command allows you to test a model using various options to customize the testing process.

    The command performs the following steps:
    - Sets up the model tester with the specified model and options.
    - Executes tests to evaluate the model's performance, metadata check, output consistency and more.
    - Generates and outputs performance metrics.

    Examples
    --------
    .. code-block:: console

        With default settings:
        $ ersilia test my_model -d /path/to/model

        With deep testing level and inspect:
        $ ersilia test my_model -d /path/to/model --level deep --inspect --remote
    """

    @ersilia_cli.command(
        short_help="Test a model",
        help="""
        Test a local models that are under development as well as on deployment and obtain a detailed report on its expected behavior and performance
        """,
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "-l",
        "--level",
        "level",
        help="Level of testing, None: for default, deep: for deep testing",
        required=False,
        default=None,
        type=click.STRING,
    )
    @click.option(
        "-d",
        "--dir",
        "dir",
        help="Model directory",
        required=False,
        default=None,
        type=click.STRING,
    )
    @click.option(
        "--inspect",
        help="Inspect the model: More on the docs",
        is_flag=True,
        default=False,
    )
    @click.option(
        "--remote",
        help="Test the model from remote git repository",
        is_flag=True,
        default=False,
    )
    @click.option(
        "--remove",
        help="Remove the model directory after testing",
        is_flag=True,
        default=False,
    )
    def test(model, level, dir, inspect, remote, remove):
        mt = ModelTester(
            model_id=model,
            level=level,
            dir=dir,
            inspect=inspect,
            remote=remote,
            remove=remove,
        )
        echo("Setting up model tester...")
        mt.setup()
        echo("Testing model...")
        mt.run(output_file=None)
