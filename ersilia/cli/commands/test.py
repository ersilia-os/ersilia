import click
from ...cli import echo
from . import ersilia_cli
from ...publish.test import ModelTester

def test_cmd():

    # Example usage: ersilia test {MODEL}
    @ersilia_cli.command(
        short_help="Test a model",
        help="Test a model and obtain performance metrics",
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "-e", 
        "--env", 
        "env", 
        required=False, 
        default=None, 
        type=click.STRING
    )
    @click.option(
        "-t", 
        "--type", 
        "type", 
        required=False, 
        default=None, 
        type=click.STRING
    )
    @click.option(
        "-l", 
        "--level", 
        "level", 
        required=False, 
        default=None, 
        type=click.STRING
    )
    @click.option(
        "-d", 
        "--dir", 
        "dir", 
        required=False, 
        default=None, 
        type=click.STRING
    )
    @click.option(
        "-o", 
        "--output", 
        "output", 
        required=False, 
        default=None, 
        type=click.STRING
    )
    def test(
        model, 
        env, 
        type, 
        level,
        dir, 
        output
      ):
        mt = ModelTester(
            model_id=model,
            env=env, 
            type=type, 
            level=level, 
            dir=dir
        )
        mt.run(output)  
