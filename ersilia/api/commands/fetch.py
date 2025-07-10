import asyncio

import nest_asyncio

from ... import ModelBase, logger
from ...hub.fetch.fetch import ModelFetcher

nest_asyncio.apply()


def _fetch(mf, model_id):
    """
    Fetches an Ersilia model to run it locally.

    This command allows users to fetch a specified model from the model hub (dockerhub, repo, s3 etc...).

    Returns
    -------
    Function: The fetch command function to be used by the API. 
    Str: Confirmation message on success or warning message on failure.
    
    Raises
    -------
    RuntimeError: If both BentoML and FastAPI are used togehter. 

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
    with_bentoml,
    with_fastapi,
    verbose,
):
    if verbose:
        logger.set_verbosity(1)
    else:
        logger.set_verbosity(0)

    if with_bentoml and with_fastapi:
        raise Exception("Cannot use both BentoML and FastAPI")
    if from_dir is not None:
        mdl = ModelBase(repo_path=from_dir)
    else:
        mdl = ModelBase(model_id_or_slug=model)
    model_id = mdl.model_id
    print(f"\033[34m⬇️ Fetching model {model_id}: {mdl.slug}\033[0m")
    mf = ModelFetcher(
        repo_path=from_dir,
        overwrite=overwrite,
        force_from_github=from_github,
        force_from_s3=from_s3,
        force_from_dockerhub=from_dockerhub,
        img_version=version,
        force_from_hosted=from_hosted,
        force_with_bentoml=with_bentoml,
        force_with_fastapi=with_fastapi,
        hosted_url=hosted_url,
        local_dir=from_dir,
    )
    fetch_result = _fetch(mf, model_id)

    if fetch_result.fetch_success:
        print(f"\033[32m👍Model {model_id} fetched successfully!\033[0m")
    else:
        print(
            f"\033[31m👎 Model {model_id} failed to fetch! {fetch_result.reason}\033[0m"
        )
