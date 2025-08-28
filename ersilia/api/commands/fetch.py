import asyncio

import nest_asyncio

from ... import ModelBase, logger
from ...hub.fetch.fetch import ModelFetcher
from ..echo import echo

nest_asyncio.apply()


def _fetch(mf, model_id):
    """
    Fetches an Ersilia model to run it locally.

    This command allows users to fetch a specified model from the model hub (dockerhub, repo, s3 etc...).

    Returns
    -------
    Function: The fetch command function to be used by the API.
    Str: Confirmation message on success or warning message on failure.

    # Raises
    # -------
    # RuntimeError: If both BentoML and FastAPI are used together.

    """
    res = asyncio.run(mf.fetch(model_id))
    return res


def fetch(
    model,
    overwrite,
    from_dir,
    from_github,
    from_dockerhub,
    version,
    from_s3,
    from_hosted,
    hosted_url,
    verbose_flag,
):
    if verbose_flag:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)

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
        overwrite=overwrite,
        force_from_github=from_github,
        force_from_s3=from_s3,
        force_from_dockerhub=from_dockerhub,
        img_version=version,
        force_from_hosted=from_hosted,
        hosted_url=hosted_url,
        local_dir=from_dir,
    )
    fetch_result = _fetch(mf, model_id)

    if fetch_result.fetch_success:
        echo(
            ":thumbs_up: Model {0} fetched successfully!".format(model_id),
            fg="green",
        )
    elif (
        fetch_result.reason
        == "Model already exists on your system. If you want to fetch it again, please delete the existing model first."
    ):
        echo(
            ":thumbs_up: Model {0} is already fetched successfully!".format(model_id),
            fg="green",
        )
    else:
        echo(
            f":thumbs_down: Model {model_id} failed to fetch! {fetch_result.reason}",
            fg="red",
        )
        echo(fetch_result.reason)
