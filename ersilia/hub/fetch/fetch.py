"""Fetch Model from the Ersilia Model Hub."""
import os
import json
import importlib

from .lazy_fetchers.dockerhub import ModelDockerHubFetcher
from .lazy_fetchers.hosted import ModelHostedFetcher
from .register.standard_example import ModelStandardExample
from ... import ErsiliaBase

from . import STATUS_FILE, DONE_TAG


class ModelFetcher(ErsiliaBase):
    def __init__(
        self,
        config_json=None,
        credentials_json=None,
        overwrite=None,
        repo_path=None,
        mode=None,
        pip=False,
        dockerize=False,
        force_from_github=False,
        force_from_s3=False,
        force_from_dockerhub=False,
        force_from_hosted=False,
        force_with_bentoml=False,
        hosted_url=None,
    ):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.overwrite = overwrite
        self.mode = mode
        self.do_pip = pip
        self.repo_path = repo_path
        if self.mode == "docker":
            self.logger.debug("When packing mode is docker, dockerization is mandatory")
            dockerize = True
        self.do_docker = dockerize
        self.model_dockerhub_fetcher = ModelDockerHubFetcher(
            overwrite=self.overwrite, config_json=self.config_json
        )
        self.is_docker_installed = self.model_dockerhub_fetcher.is_docker_installed()
        self.is_docker_active = self.model_dockerhub_fetcher.is_docker_active()
        self.model_hosted_fetcher = ModelHostedFetcher(
            url=hosted_url, config_json=self.config_json
        )
        self.force_from_github = force_from_github
        self.force_from_s3 = force_from_s3
        self.force_from_dockerhub = force_from_dockerhub
        self.force_from_hosted = force_from_hosted
        self.force_with_bentoml = force_with_bentoml
        self.hosted_url = hosted_url

    def _fetch_not_from_dockerhub(self, model_id):
        self.model_id = model_id
        if not self.force_with_bentoml:
            self.logger.debug("Fetching using Ersilia Pack (FastAPI)")
            fetch = importlib.import_module("ersilia.hub.fetch.fetch_fastapi")
            mf = fetch.ModelFetcherFromFastAPI()
            if mf.seems_installable():
                mf.fetch()
        if not self.exists(model_id):
            self.logger.debug("Fetching using BentoML")
            fetch = importlib.import_module("ersilia.hub.fetch.fetch_bentoml")
            mf = fetch.ModelFetcherFromBentoML()
            if mf.seems_installable():
                mf.fetch()

    def _fetch_from_dockerhub(self, model_id):
        self.logger.debug("Fetching from DockerHub")
        self.model_dockerhub_fetcher.fetch(model_id=model_id)

    def _fetch_from_hosted(self, model_id):
        self.logger.debug("Fetching from hosted")
        self.model_hosted_fetcher.fetch(model_id=model_id)
        self.logger.debug("Fetching from hosted done")

    def _decide_if_use_dockerhub(self, model_id):
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
            self.logger.debug("Docker is not installed in your local")
            return False
        if not self.is_docker_active:
            self.logger.debug("Docker is not active in your local")
            return False
        if not self.model_dockerhub_fetcher.is_available(model_id=model_id):
            self.logger.debug("Docker image of this model doesn't seem to be available")
            return False
        return True

    def _decide_if_use_hosted(self, model_id):
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

    def exists(self, model_id):
        status_file = os.path.join(self._model_path(model_id), STATUS_FILE)
        if not os.path.exists(status_file):
            return False
        with open(status_file, "r") as f:
            status = json.load(f)
        if status["done"]:
            return True
        else:
            return False

    def _standard_csv_example(self, model_id):
        ms = ModelStandardExample(model_id=model_id, config_json=self.config_json)
        ms.run()

    def _fetch(self, model_id):
        self.logger.debug("Starting fetching procedure")
        do_hosted = self._decide_if_use_hosted(model_id=model_id)
        if do_hosted:
            self._fetch_from_hosted(model_id=model_id)
            return
        do_dockerhub = self._decide_if_use_dockerhub(model_id=model_id)
        if do_dockerhub:
            self._fetch_from_dockerhub(model_id=model_id)
            return
        if self.overwrite is None:
            self.overwrite = True
        self._fetch_not_from_dockerhub(model_id=model_id)

    def fetch(self, model_id):
        self._fetch(model_id)
        self._standard_csv_example(model_id)
