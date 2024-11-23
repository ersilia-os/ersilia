import click
import asyncio
import nest_asyncio
from . import ersilia_cli
from .. import echo
from ...hub.fetch.fetch import ModelFetcher
from ... import ModelBase

nest_asyncio.apply()

def fetch_cmd():
    """Create fetch commmand"""

    def _fetch(mf, model_id):
        res = asyncio.run(mf.fetch(model_id))
        return res

    # Example usage: ersilia fetch {MODEL}
    @ersilia_cli.command(
        short_help="Fetch model from Ersilia Model Hub",
        help="Fetch model from EOS repository. Model files are downloaded from GitHub and model data are "
        "downloaded from a file storage system such as the Open Science Framework. Model is downloaded to "
        "an EOS folder, then packed to a BentoML bundle",
    )
    @click.argument("model", type=click.STRING)
    @click.option("--mode", "-m", default=None, type=click.STRING)
    @click.option("--dockerize/--not-dockerize", default=False)
    @click.option(
        "--overwrite/--reuse",
        default=True,
        help="Overwrite environment or reuse using already available environment for this model",
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
        "--from_hosted",
        is_flag=True,
        default=False,
        help="Force fetch from hosted service. This only creates a basic folder structure for the model, the model is not actually downloaded.",
    )
    @click.option(
        "--from_url",
        type=click.STRING,
        default=None,
        help="Fetch a model based on a URL. This only creates a basic folder structure for the model, the model is not actually downloaded.",
    )
    @click.option(
        "--with_bentoml",
        is_flag=True,
        default=False,
        help="Force fetch using BentoML",
    )
    @click.option(
        "--with_fastapi",
        is_flag=True,
        default=False,
        help="Force fetch using FastAPI",
    )
    def fetch(
        model,
        mode,
        dockerize,
        overwrite,
        from_dir,
        from_github,
        from_dockerhub,
        from_s3,
        from_hosted,
        from_url,
        with_bentoml,
        with_fastapi,
    ):
        if with_bentoml and with_fastapi:
            raise Exception("Cannot use both BentoML and FastAPI")

        if from_dir is not None:
            mdl = ModelBase(repo_path=from_dir)
        else:
            mdl = ModelBase(model_id_or_slug=model)
        model_id = mdl.model_id
        echo(
            ":down_arrow:  Fetching model {0}: {1}".format(model_id, mdl.slug),
            fg="blue",
        )
        mf = ModelFetcher(
            repo_path=from_dir,
            mode=mode,
            dockerize=dockerize,
            overwrite=overwrite,
            force_from_github=from_github,
            force_from_s3=from_s3,
            force_from_dockerhub=from_dockerhub,
            force_from_hosted=from_hosted,
            force_with_bentoml=with_bentoml,
            force_with_fastapi=with_fastapi,
            hosted_url=from_url,
            local_dir=from_dir,
        )
        is_fetched = _fetch(mf, model_id)

        if is_fetched:
            echo(":thumbs_up: Model {0} fetched successfully!".format(model_id), fg="green")
        else:
            echo(":thumbs_down: Model {0} failed to fetch!".format(model_id), fg="red")
