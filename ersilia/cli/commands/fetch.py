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
    @click.option(
        "--overwrite/--reuse",
        default=True,
        help="Overwrite environment or reuse using already available environment for this model",
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
    def fetch(
        model,
        repo_path,
        mode,
        dockerize,
        overwrite,
        from_github,
        from_dockerhub,
        from_s3,
    ):
        if repo_path is not None:
            mdl = ModelBase(repo_path=repo_path)
        else:
            mdl = ModelBase(model_id_or_slug=model)
        model_id = mdl.model_id
        echo(
            ":down_arrow:  Fetching model {0}: {1}".format(model_id, mdl.slug),
            fg="blue",
        )
        mf = ModelFetcher(
            repo_path=repo_path,
            mode=mode,
            dockerize=dockerize,
            overwrite=overwrite,
            force_from_github=from_github,
            force_from_s3=from_s3,
            force_from_dockerhub=from_dockerhub,
        )
        _fetch(mf, model_id)
        echo(":thumbs_up: Model {0} fetched successfully!".format(model_id), fg="green")
