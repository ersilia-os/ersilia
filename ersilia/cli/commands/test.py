import click
from ...cli import echo
from . import ersilia_cli
from ...publish.test import ModelTester

def test_cmd():

# Example usage: ersilia test {MODEL} -d local_dir -l deep --remote
    @ersilia_cli.command(
        short_help="Test a model",
        help="Test a model and obtain performance metrics",
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "-l", 
        "--level", 
        "level", 
        help="Level of testing, None: for default, deep: for deep testing",
        required=False, 
        default=None, 
        type=click.STRING
    )
    @click.option(
        "-d", 
        "--dir", 
        "dir", 
        help="Model directory",
        required=False, 
        default=None, 
        type=click.STRING
    )
    @click.option(
        "--inspect", 
        help="Inspect the model: More on the docs", 
        is_flag=True, 
        default=False
    )
    @click.option(
        "--remote",
        help="Test the model from remote git repository", 
        is_flag=True, 
        default=False
    )
    @click.option(
        "--remove",
        help="Remove the model after testing", 
        is_flag=True, 
        default=False
    )
    def test(
        model, 
        level,
        dir, 
        inspect,
        remote,
        remove
      ):
        mt = ModelTester(
            model_id=model,
            level=level, 
            dir=dir,
            inspect=inspect,
            remote=remote,
            remove=remove
        )
        echo("Setting up model tester...")
        mt.setup()
        echo("Testing model...")
        mt.run(output_file=None)  
