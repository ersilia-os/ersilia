import click

from ...cli import echo
from ...publish.test import ModelTester
from . import ersilia_cli


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

        With basic testing:/
        $ ersilia test eosxxxx --from_dir /path/to/model

        With different sources to fetch the model:
        $ ersilia test eosxxxx --from_github/--from_dockerhub/--from_s3

        With different levels of testing:
        $ ersilia test eosxxxx --shallow --from_github/--from_dockerhub/--from_s3
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
        help="Level of testing, None: for default, deep: for deep testing, shallow: for shallow testing",
        required=False,
        default=None,
        type=click.STRING,
    )
    @click.option(
        "--from_dir",
        default=None,
        type=click.STRING,
        help="Local path where the model is stored",
    )
    @click.option(
        "--from_github",
        is_flag=True,
        default=False,
        help="Fetch fetch directly from GitHub",
    )
    @click.option(
        "--from_dockerhub",
        is_flag=True,
        default=False,
        help="Force fetch from DockerHub",
    )
    @click.option(
        "--from_s3", is_flag=True, default=False, help="Force fetch from AWS S3 bucket"
    )
    @click.option(
        "--version",
        default=None,
        type=click.STRING,
        help="Version of the model to fetch, when fetching a model from DockerHub",
    )
    @click.option(
        "--shallow",
        is_flag=True,
        default=False,
        help="This flag is used to check shallow checks (such as container size, output consistency..)",
    )
    @click.option(
        "--deep",
        is_flag=True,
        default=False,
        help="This flag is used to check deep checks (such as computational performance checks)",
    )
    @click.option(
        "--as_json",
        is_flag=True,
        default=False,
        help="This flag is used to save the report as json file)",
    )
    def test(
        model,
        level,
        from_dir,
        from_github,
        from_dockerhub,
        from_s3,
        version,
        shallow,
        deep,
        as_json,
    ):
        mt = ModelTester(
            model,
            level,
            from_dir,
            from_github,
            from_dockerhub,
            from_s3,
            version,
            shallow,
            deep,
            as_json,
        )
        echo(f"Model testing started for: {model}")
        mt.run()
