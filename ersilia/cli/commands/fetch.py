import click
from . import ersilia_cli
from .. import echo
from ...hub.fetch.fetch import ModelFetcher
from ... import ModelBase


def fetch_cmd():
    """Create fetch commmand"""

    def _fetch(mf, model_id):
        mf.fetch(model_id)

    # Example usage: ersilia fetch {MODEL}
    @ersilia_cli.command(
        short_help="Fetch model from Ersilia Model Hub",
        help="Fetch model from EOS repository. Model files are downloaded from GitHub and model data are "
        "downloaded from a file storage system such as the Open Science Framework. Model is downloaded to "
        "an EOS folder, then packed to a BentoML bundle",
    )
    @click.argument("model", type=click.STRING)
    @click.option("--repo_path", "-r", default=None, type=click.STRING)
    @click.option("--mode", "-m", default=None, type=click.STRING)
    @click.option("--dockerize/--not-dockerize", default=False)
    def fetch(model, repo_path, mode, dockerize):
        mdl = ModelBase(model)
        model_id = mdl.model_id
        echo(
            ":down_arrow:  Fetching model {0}: {1}".format(model_id, mdl.slug),
            fg="blue",
        )
        mf = ModelFetcher(repo_path=repo_path, mode=mode, dockerize=dockerize)
        _fetch(mf, model_id)
        echo(":thumbs_up: Model {0} fetched successfully!".format(model_id), fg="green")
