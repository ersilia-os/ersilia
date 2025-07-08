
import asyncio
import nest_asyncio
from ... import ModelBase
from ...hub.fetch.fetch import ModelFetcher
nest_asyncio.apply()
 
def _fetch(mf, model_id):
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
):
    if with_bentoml and with_fastapi:
            raise Exception("Cannot use both BentoML and FastAPI")
    if from_dir is not None:
        mdl = ModelBase(repo_path=from_dir)
    else:
        mdl = ModelBase(model_id_or_slug=model)
    model_id = mdl.model_id
    print(
        ":down_arrow:  Fetching model {0}: {1}".format(model_id, mdl.slug),
    )
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
        print(
            ":thumbs_up: Model {0} fetched successfully!".format(model_id)
        )
    else:
        print(
            f":thumbs_down: Model {model_id} failed to fetch! {fetch_result.reason}"
        )