import click

from ...cli import echo
from ...publish.test.test import ModelTester
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
        help="This flag is used to test shallow checks (fetch and run models and output consistency.)",
    )
    @click.option(
        "--deep",
        is_flag=True,
        default=False,
        help="This flag is used to test deep checks (computational performance)",
    )
    @click.option(
        "--surface",
        is_flag=True,
        default=False,
        help="This flag is used to test surface checks (simple model fetch and run)",
    )
    @click.option(
        "--inspect",
        is_flag=True,
        default=False,
        help="This flag is used to check inspect checks (file structure and metadata)",
    )
    @click.option(
        "--report_path",
        default=None,
        type=click.STRING,
        help="This flag is used to specify a path for report json file.)",
    )
    @click.option(
        "--clean",
        is_flag=True,
        default=False,
        help="This flag is used to clean out the temp folder after testing",
    )
    def test(
        model,
        from_dir,
        from_github,
        from_dockerhub,
        from_s3,
        version,
        shallow,
        deep,
        surface,
        inspect,
        report_path,
        clean,
    ):
        mt = ModelTester(
            model,
            from_dir,
            from_github,
            from_dockerhub,
            from_s3,
            version,
            shallow,
            deep,
            surface,
            inspect,
            report_path,
            clean,
        )
        echo(f"Model testing started for: {model}")
        mt.run()

    return test
