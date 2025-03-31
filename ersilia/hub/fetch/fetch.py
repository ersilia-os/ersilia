import importlib
import json
import os
from collections import namedtuple

from ... import ErsiliaBase
from ...db.hubdata.interfaces import JsonModelsInterface
from ...default import (
    MODEL_SOURCE_FILE,
    PACK_METHOD_BENTOML,
    PACK_METHOD_FASTAPI,
)
from ...hub.delete.delete import ModelFullDeleter
from ...hub.fetch.actions.template_resolver import TemplateResolver
from ...setup.requirements import check_bentoml
from ...tools.bentoml.exceptions import BentoMLException
from ...utils.exceptions_utils.fetch_exceptions import (
    NotInstallableWithBentoML,
    NotInstallableWithFastAPI,
    StandardModelExampleError,
    WithToolFetchingNotWorking,
)
from ...utils.exceptions_utils.throw_ersilia_exception import throw_ersilia_exception
from ...utils.terminal import yes_no_input
from . import STATUS_FILE
from .lazy_fetchers.dockerhub import ModelDockerHubFetcher
from .lazy_fetchers.hosted import ModelHostedFetcher
from .register.standard_example import ModelStandardExample

FetchResult = namedtuple("FetchResult", ["fetch_success", "reason"])


class ModelFetcher(ErsiliaBase):
    """
    ModelFetcher is responsible for fetching models from various sources.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings for the fetcher.
    credentials_json : dict, optional
        Credentials for accessing the model hub.
    overwrite : bool, optional
        Whether to overwrite existing files.
    repo_path : str, optional
        Path to the repository.
    mode : str, optional
        Mode of operation, e.g., 'docker'.
    pip : bool, optional
        Whether to use pip for installation.
    dockerize : bool, optional
        Whether to dockerize the model.
    force_from_github : bool, optional
        Whether to force fetching from GitHub.
    force_from_s3 : bool, optional
        Whether to force fetching from S3.
    force_from_dockerhub : bool, optional
        Whether to force fetching from DockerHub.
    img_version : str, optional
        Version of the model image.
    force_from_hosted : bool, optional
        Whether to force fetching from hosted services.
    force_with_bentoml : bool, optional
        Whether to force fetching with BentoML.
    force_with_fastapi : bool, optional
        Whether to force fetching with FastAPI.
    hosted_url : str, optional
        URL for hosted model.
    local_dir : str, optional
        Local directory for the model.

    Examples
    --------
    .. code-block:: python

        fetcher = ModelFetcher(config_json=config)
        await fetcher.fetch(model_id="eosxxxx")
    """

    def __init__(
        self,
        config_json: dict = None,
        credentials_json: dict = None,
        overwrite: bool = None,
        repo_path: str = None,
        mode: str = None,
        pip: bool = False,
        dockerize: bool = False,
        force_from_github: bool = False,
        force_from_s3: bool = False,
        force_from_dockerhub: bool = False,
        img_version: str = None,
        force_from_hosted: bool = False,
        force_with_bentoml: bool = False,
        force_with_fastapi: bool = False,
        hosted_url: str = None,
        local_dir: str = None,
    ):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )

        self.ji = JsonModelsInterface(config_json=self.config_json)
        self.overwrite = overwrite
        self.mode = mode
        self.do_pip = pip
        self.repo_path = repo_path
        if self.mode == "docker":
            self.logger.debug("When packing mode is docker, dockerization is mandatory")
            dockerize = True
        self.do_docker = dockerize
        self.model_dockerhub_fetcher = ModelDockerHubFetcher(
            overwrite=self.overwrite,
            config_json=self.config_json,
            img_tag=img_version,
            force_with_bentoml=force_with_bentoml,
            force_with_fastapi=force_with_fastapi,
        )
        self.is_docker_installed = self.model_dockerhub_fetcher.is_docker_installed()
        self.is_docker_active = self.model_dockerhub_fetcher.is_docker_active()
        self.model_hosted_fetcher = ModelHostedFetcher(
            url=hosted_url, config_json=self.config_json
        )
        self.check_bentoml = check_bentoml
        self.can_use_docker = self.is_docker_installed and self.is_docker_active
        self.force_from_github = force_from_github
        self.force_from_s3 = force_from_s3
        self.force_from_dockerhub = force_from_dockerhub
        self.force_from_hosted = force_from_hosted
        self.force_with_bentoml = force_with_bentoml
        self.force_with_fastapi = force_with_fastapi
        self.hosted_url = hosted_url
        self.local_dir = local_dir

        self.logger.debug("Getting model source")
        sources = {
            self.force_from_github: "GitHub",
            self.force_from_s3: "Amazon S3",
            self.force_from_dockerhub: "DockerHub",
            self.force_from_hosted: "Hosted services",
            self.local_dir is not None: "Local Repository",
        }

        self.model_source = next(
            (source for condition, source in sources.items() if condition), "DockerHub"
        )
        self.logger.debug("Model getting fetched from {0}".format(self.model_source))

    @throw_ersilia_exception()
    def _decide_fetcher(self, model_id: str) -> str:
        tr = TemplateResolver(model_id=model_id, repo_path=self.repo_path)
        if tr.is_fastapi():
            return PACK_METHOD_FASTAPI
        elif tr.is_bentoml():
            return PACK_METHOD_BENTOML
        else:
            raise Exception("No fetcher available")

    @throw_ersilia_exception()
    def _fetch_from_fastapi(self):
        self.logger.debug("Fetching using Ersilia Pack (FastAPI)")
        fetch = importlib.import_module("ersilia.hub.fetch.fetch_fastapi")
        mf = fetch.ModelFetcherFromFastAPI(
            config_json=self.config_json,
            credentials_json=self.credentials_json,
            overwrite=self.overwrite,
            repo_path=self.repo_path,
            mode=self.mode,
            force_from_github=self.force_from_github,
            force_from_s3=self.force_from_s3,
        )
        if mf.seems_installable(model_id=self.model_id):
            mf.fetch(model_id=self.model_id)
        else:
            self.logger.debug("Not installable with FastAPI")
            raise NotInstallableWithFastAPI(model_id=self.model_id)

    @throw_ersilia_exception()
    def _fetch_from_bentoml(self):
        self.logger.debug("Fetching using BentoML")
        self.check_bentoml()

        fetch = importlib.import_module("ersilia.hub.fetch.fetch_bentoml")
        mf = fetch.ModelFetcherFromBentoML(
            config_json=self.config_json,
            credentials_json=self.credentials_json,
            overwrite=self.overwrite,
            repo_path=self.repo_path,
            mode=self.mode,
            pip=self.do_pip,
            dockerize=self.do_docker,
            force_from_github=self.force_from_github,
            force_from_s3=self.force_from_s3,
        )

        # Check if the model can be installed with BentoML
        if mf.seems_installable(model_id=self.model_id):
            mf.fetch(model_id=self.model_id)
        else:
            raise NotInstallableWithBentoML(model_id=self.model_id)

    @throw_ersilia_exception()
    def _fetch_not_from_dockerhub(self, model_id: str):
        self.model_id = model_id
        is_fetched = False

        tr = TemplateResolver(model_id=model_id, repo_path=self.repo_path)
        if self.force_with_fastapi:
            if not tr.is_fastapi():
                raise WithToolFetchingNotWorking(tool="fastapi")
            self._fetch_from_fastapi()
            is_fetched = True
        if self.force_with_bentoml:
            if not tr.is_bentoml():
                raise WithToolFetchingNotWorking(tool="bentoml")
            self._fetch_from_bentoml()
            is_fetched = True

        if is_fetched:
            return
        else:
            self.logger.debug("Deciding fetcher (BentoML or FastAPI)")
            fetcher_type = self._decide_fetcher(model_id)
            if fetcher_type == PACK_METHOD_FASTAPI:
                self._fetch_from_fastapi()
            if fetcher_type == PACK_METHOD_BENTOML:
                self._fetch_from_bentoml()

        self.logger.debug("Model already exists in your local, skipping fetching")

    def _standard_csv_example(self, model_id: str):
        ms = ModelStandardExample(model_id=model_id, config_json=self.config_json)
        ms.run()

    async def _fetch_from_dockerhub(self, model_id: str):
        self.logger.debug("Fetching from DockerHub")
        await self.model_dockerhub_fetcher.fetch(model_id=model_id)

    def _fetch_from_hosted(self, model_id: str):
        self.logger.debug("Fetching from hosted")
        self.model_hosted_fetcher.fetch(model_id=model_id)
        self.logger.debug("Fetching from hosted done")

    def _decide_if_use_dockerhub(self, model_id: str) -> bool:
        if self.repo_path is not None:
            return False
        if self.force_from_dockerhub:
            return True
        if self.force_from_s3:
            return False
        if self.force_from_github:
            return False
        if self.force_from_hosted:
            return False
        if not self.is_docker_installed:
            self.logger.debug("Docker Engine is not installed on your system.")
            return False
        if not self.is_docker_active:
            self.logger.info("Docker is not active in your local")
            return False
        if not self.ji.identifier_exists(model_id=model_id):
            self.logger.info("Docker image of this model doesn't seem to be available")
            return False
        return True

    def _decide_if_use_hosted(self, model_id: str) -> bool:
        if self.repo_path is not None:
            return False
        if self.force_from_dockerhub:
            return False
        if self.force_from_github:
            return False
        if self.force_from_s3:
            return False
        if not self.model_hosted_fetcher.is_available(model_id=model_id):
            self.logger.debug("There is no hosted URL available for this model")
            return False
        if self.force_from_hosted:
            return True
        return False

    def exists(self, model_id: str) -> bool:
        """
        Check if the model exists locally.

        Parameters
        ----------
        model_id : str
            The ID of the model to be checked.

        Returns
        -------
        bool
            True if the model exists locally, False otherwise.
        """
        status_file = os.path.join(self._model_path(model_id), STATUS_FILE)
        if not os.path.exists(status_file):
            return False
        with open(status_file, "r") as f:
            status = json.load(f)
        if status["done"]:
            return True
        else:
            return False

    async def _fetch(self, model_id: str) -> FetchResult:
        if not self.exists(model_id):
            self.logger.info("Model doesn't exist on your system, fetching it now.")
            self.logger.debug("Starting fetching procedure")
            do_dockerhub = self._decide_if_use_dockerhub(model_id=model_id)
            if do_dockerhub:
                self.logger.debug("Decided to fetch from DockerHub")
                if not self.can_use_docker:
                    return FetchResult(
                        fetch_success=False,
                        reason="Docker is not installed or active on your system.",
                    )
                await self._fetch_from_dockerhub(model_id=model_id)
                return FetchResult(
                    fetch_success=True, reason="Model fetched successfully"
                )
            do_hosted = self._decide_if_use_hosted(model_id=model_id)
            if do_hosted:
                self.logger.debug("Fetching from hosted")
                self._fetch_from_hosted(model_id=model_id)
                return FetchResult(
                    fetch_success=True, reason="Model fetched successfully"
                )
            if self.overwrite is None:
                self.logger.debug("Overwriting")
                self.overwrite = True
            self.logger.debug("Fetching in your system, not from DockerHub")
            self._fetch_not_from_dockerhub(model_id=model_id)
            return FetchResult(fetch_success=True, reason="Model fetched successfully")
        else:
            self.logger.info(
                "Model already exists on your system. If you want to fetch it again, please delete it first."
            )
            return FetchResult(
                fetch_success=False,
                reason="Model already exists on your system. If you want to fetch it again, please delete the existing model first.",
            )

    async def fetch(self, model_id: str) -> bool:
        """
        Fetch a model with the given eos identifier.

        Parameters
        ----------
        model_id : str
            The eos identifier of the model.

        Returns
        -------
        bool
            True if the model was fetched successfully, False otherwise.

        Examples
        --------
        .. code-block:: python

            fetcher = ModelFetcher(config_json=config)
            success = await fetcher.fetch(model_id="eosxxxx")
        """
        try:
            fr = await self._fetch(model_id)
            if not fr.fetch_success:
                return fr

            self._standard_csv_example(model_id)
            self.logger.debug("Writing model source to file")
            model_source_file = os.path.join(
                self._model_path(model_id), MODEL_SOURCE_FILE
            )
            try:
                os.makedirs(self._model_path(model_id), exist_ok=True)
            except OSError as error:
                self.logger.error(f"Error during folder creation: {error}")
            with open(model_source_file, "w") as f:
                f.write(self.model_source)

            return FetchResult(fetch_success=True, reason="Model fetched successfully")

        except (StandardModelExampleError, BentoMLException) as err:
            self.logger.debug(f"{type(err).__name__} occurred: {str(err)}")
            do_delete = yes_no_input(
                "Do you want to delete the model artifacts? [Y/n]",
                default_answer="Y",
            )
            if do_delete:
                md = ModelFullDeleter(overwrite=False)
                md.delete(model_id)
                self.logger.info(
                    f"✅ Model '{model_id}' artifacts have been successfully deleted."
                )

            reason = (
                str(err) if str(err) else "An unknown error occurred during fetching."
            )
            return FetchResult(fetch_success=False, reason=reason)
