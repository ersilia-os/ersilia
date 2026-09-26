import asyncio
import sys

import rich_click as click

from .. import echo
from . import ersilia_cli


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
        default=False,
        help="Fetch from DockerHub (the default when no other source is given).",
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
        import nest_asyncio

        from ... import ModelBase
        from ...hub.fetch.fetch import ALREADY_FETCHED, ModelFetcher

        nest_asyncio.apply()

        from ..run_checks import fail

        chosen = [
            name
            for name, on in (
                ("--from_dir", from_dir),
                ("--from_github", from_github),
                ("--from_s3", from_s3),
                ("--from_dockerhub", from_dockerhub),
                ("--from_hosted", from_hosted),
            )
            if on
        ]
        if len(chosen) > 1:
            raise click.UsageError(
                "Choose only one source; got {0}.".format(", ".join(chosen))
            )
        from_dockerhub = not chosen or from_dockerhub
        if version is not None and not from_dockerhub:
            echo("--version only applies to DockerHub, so it is ignored.", fg="yellow")

        if from_dir is not None:
            import os

            if not os.path.isdir(os.path.expanduser(from_dir)):
                fail(f"The folder {from_dir} does not exist.")
            mdl = ModelBase(repo_path=from_dir)
            from ...utils.paths import get_metadata_from_base_dir

            try:
                folder_id = get_metadata_from_base_dir(from_dir).get("Identifier")
            except Exception:
                folder_id = None
            folder_id = folder_id or mdl.model_id
            requested = model.strip().lower()
            if folder_id and requested not in (folder_id, mdl.slug):
                fail(
                    f"The folder {from_dir} contains model {folder_id}, not {model}.",
                    f"Run 'ersilia fetch {folder_id} --from_dir {from_dir}'.",
                )
        else:
            mdl = ModelBase(model_id_or_slug=model)
        model_id = mdl.model_id

        mf = ModelFetcher(
            repo_path=from_dir,
            force_from_github=from_github,
            force_from_s3=from_s3,
            force_from_dockerhub=from_dockerhub,
            img_version=version,
            force_from_hosted=from_hosted is not None,
            hosted_url=from_hosted,
            local_dir=from_dir,
            slug=mdl.slug,
        )
        fetch_result = _fetch(mf, model_id)

        if fetch_result.fetch_success:
            if fetch_result.reason == "Model fetched successfully":
                echo(f"Model {model_id} fetched.", fg="green")
        elif fetch_result.reason == ALREADY_FETCHED:
            echo(f"Model {model_id} is already fetched.", fg="yellow")
            echo(
                f"To fetch it again, delete it first with 'ersilia delete {model_id}'."
            )
        else:
            echo(f"Model {model_id} could not be fetched.", fg="red")
            echo(fetch_result.reason)
            sys.exit(1)

    return fetch
