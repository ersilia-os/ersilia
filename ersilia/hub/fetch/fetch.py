"""Fetch model from the Ersilia Model Hub"""

import json
import os
from timeit import default_timer as timer
from datetime import timedelta

from .actions.setup import SetupChecker
from .actions.prepare import ModelPreparer
from .actions.get import ModelGetter
from .actions.lake import LakeGetter
from .actions.pack import ModelPacker
from .actions.toolize import ModelToolizer
from .actions.content import CardGetter
from .actions.check import ModelChecker
from .actions.sniff import ModelSniffer
from .actions.inform import ModelInformer
from .register.register import ModelRegisterer
from .lazy_fetchers.dockerhub import ModelDockerHubFetcher
from .lazy_fetchers.hosted import ModelHostedFetcher

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
        self.model_hosted_fetcher = ModelHostedFetcher(
            url=hosted_url, config_json=self.config_json
        )
        self.force_from_github = force_from_github
        self.force_from_s3 = force_from_s3
        self.force_from_dockerhub = force_from_dockerhub
        self.force_from_hosted = force_from_hosted
        self.hosted_url = hosted_url

    def _setup_check(self):
        sc = SetupChecker(model_id=self.model_id, config_json=self.config_json)
        sc.check()

    def _prepare(self):
        mp = ModelPreparer(
            model_id=self.model_id,
            overwrite=self.overwrite,
            config_json=self.config_json,
        )
        mp.prepare()

    def _get(self):
        mg = ModelGetter(
            model_id=self.model_id,
            repo_path=self.repo_path,
            config_json=self.config_json,
            force_from_gihtub=self.force_from_github,
            force_from_s3=self.force_from_s3,
        )
        mg.get()

    def _lake(self):
        ml = LakeGetter(model_id=self.model_id, config_json=self.config_json)
        ml.get()

    def _pack(self):
        mp = ModelPacker(
            model_id=self.model_id, mode=self.mode, config_json=self.config_json
        )
        mp.pack()

    def _toolize(self):
        mt = ModelToolizer(model_id=self.model_id, config_json=self.config_json)
        mt.toolize(do_pip=self.do_pip, do_docker=self.do_docker)

    def _content(self):
        cg = CardGetter(self.model_id, self.config_json)
        cg.get()

    def _check(self):
        mc = ModelChecker(self.model_id, self.config_json)
        mc.check()

    def _sniff(self):
        sn = ModelSniffer(self.model_id, self.config_json)
        sn.sniff()

    def _inform(self):
        mi = ModelInformer(self.model_id, self.config_json)
        mi.inform()

    def _success(self):
        done = {DONE_TAG: True}
        status_file = os.path.join(self._dest_dir, self.model_id, STATUS_FILE)
        with open(status_file, "w") as f:
            json.dump(done, f, indent=4)
        mr = ModelRegisterer(self.model_id, config_json=self.config_json)
        mr.register(is_from_dockerhub=False)

    def _fetch_not_from_dockerhub(self, model_id):
        start = timer()
        self.model_id = model_id
        self._setup_check()
        self._prepare()
        self._get()
        self._pack()
        self._toolize()
        self._content()
        self._check()
        self._sniff()
        self._inform()
        self._success()
        end = timer()
        elapsed_time = timedelta(seconds=end - start)
        self.logger.debug(
            "Fetching {0} done in time: {1}s".format(model_id, abs(elapsed_time))
        )
        self.logger.info(
            "Fetching {0} done successfully: {1}".format(model_id, elapsed_time)
        )

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

    def fetch(self, model_id):
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
