import click

from . import ersilia_cli
from .. import echo
from ...hub.fetch import ModelFetcher


def fetch_cmd():
    """Create fetch commmand"""
    # Example usage: ersilia fetch {MODEL_ID}
    @ersilia_cli.command(
        short_help="Fetch model from Ersilia Model Hub",
        help="Fetch model from EOS repository. Model files are downloaded from GitHub and model data are "
             "downloaded from a file storage system such as the Open Science Framework. Model is downloaded to "
             "an EOS folder, then packed to a BentoML bundle",
    )
    @click.argument("model_id", type=click.STRING)
    def fetch(model_id):
        echo("Fetching model {0}".format(model_id))
        mf = ModelFetcher()
        mf.fetch(model_id)
        echo(":thumbs_up: Model {0} fetched successfully!".format(model_id), fg="green")
