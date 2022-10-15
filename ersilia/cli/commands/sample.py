import click
from . import ersilia_cli
from .. import echo
from ...db.hubdata.samplers import InputSampler, ModelSampler
from ...core.session import Session


def sample_cmd():
    @ersilia_cli.command(
        short_help="Sample inputs and model identifiers",
        help="Sample inputs and model identifiers",
    )
    @click.option(
        "-m",
        "--model",
        is_flag=True,
        default=False,
        help="Sample a model identifier. The default is to sample inputs",
    )
    @click.option(
        "--n_samples",
        "-n",
        default=1,
        type=click.INT,
        help="Number of entities to sample",
    )
    @click.option(
        "--file_name", "-f", default=None, type=click.STRING, help="Output file name"
    )
    def sample(model, n_samples, file_name):
        session = Session(config_json=None)
        if model:
            sampler = ModelSampler(config_json=None)
            sampler.sample(n_samples=n_samples, file_name=file_name)
        else:
            model_id = session.current_model_id()
            if model_id is None:
                echo("No model is served at the moment.")
            sampler = InputSampler(model_id=model_id, config_json=None)
            sampler.sample(n_samples=n_samples, file_name=file_name)
