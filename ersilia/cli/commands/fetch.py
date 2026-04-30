import asyncio

import nest_asyncio
import rich_click as click

from ... import ModelBase
from ...hub.fetch.fetch import ModelFetcher
from .. import echo
from . import ersilia_cli

nest_asyncio.apply()


def fetch_cmd():
    """
    Fetches a specified model.

    This command allows users to fetch a specified model from the model hub (dockerhub, repo, s3 etc...).

    Returns
    -------
    function
        The fetch command function to be used by the CLI and for testing in the pytest.

    Examples
    --------
    .. code-block:: console

        Fetch a model by its ID:
        $ ersilia fetch <model_id> [auto model source decider] or ersilia fetch <model_id> --from_github/--from_dockerhub

        Fetch a model from a local directory:
        $ ersilia fetch <model_id> --from_dir <path>
    """

    def _fetch(mf, model_id):
        res = asyncio.run(mf.fetch(model_id))
        return res

    # Example usage: ersilia fetch {MODEL}
    @ersilia_cli.command(
        short_help="Fetch a model from the Ersilia Model Hub",
        help="Download a model from the Ersilia Model Hub and set it up locally. By default, models are fetched from DockerHub. Use the --from_* flags to override the source.",
    )
    @click.argument("model", type=click.STRING)
    @click.option(
        "--from_dir",
        default=None,
        type=click.STRING,
        help="Fetch from a local directory containing the model repository.",
    )
    @click.option(
        "--from_github",
        is_flag=True,
        default=False,
        help="Fetch directly from the model's GitHub repository.",
    )
    @click.option(
        "--from_s3",
        is_flag=True,
        default=False,
        help="Fetch from the Ersilia AWS S3 bucket.",
    )
    @click.option(
        "--from_dockerhub",
        is_flag=True,
        default=True,
        help="Fetch from DockerHub (default).",
    )
    @click.option(
        "--from_hosted",
        default=None,
        type=click.STRING,
        help="Connect to a remotely hosted model service by providing its URL. Only creates a local folder structure; the model is not downloaded.",
    )
    @click.option(
        "--version",
        default=None,
        type=click.STRING,
        help="Specific Docker image version to fetch from DockerHub.",
    )
    def fetch(
        model,
        from_dir,
        from_github,
        from_dockerhub,
        version,
        from_s3,
        from_hosted,
    ):
        if from_dir is not None:
            mdl = ModelBase(repo_path=from_dir)
        else:
            mdl = ModelBase(model_id_or_slug=model)
        model_id = mdl.model_id
        echo("Fetching model {0}: {1}".format(model_id, mdl.slug), fg="cyan")

        if any([from_dir, from_github, from_s3, from_hosted]):
            from_dockerhub = False

        mf = ModelFetcher(
            repo_path=from_dir,
            force_from_github=from_github,
            force_from_s3=from_s3,
            force_from_dockerhub=from_dockerhub,
            img_version=version,
            force_from_hosted=from_hosted is not None,
            hosted_url=from_hosted,
            local_dir=from_dir,
        )
        fetch_result = _fetch(mf, model_id)

        if fetch_result.fetch_success:
            echo(f"Model {model_id} fetched successfully.", fg="green")
        else:
            echo(f"Model {model_id} failed to fetch: {fetch_result.reason}", fg="red", bold=True)

    return fetch
