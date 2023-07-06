"""Fetch model from the Ersilia Model Hub"""

import json
import os
from timeit import default_timer as timer
from datetime import timedelta
import time
from tqdm import tqdm
import datetime
import shutil

from ... import ErsiliaBase
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

from ..pull.pull import ModelPuller
from ...serve.services import PulledDockerImageService
from ...setup.requirements.docker import DockerRequirement
from ...utils.docker import SimpleDocker

from . import STATUS_FILE, DONE_TAG
from ... import EOS
from ...default import (
    IS_FETCHED_FROM_DOCKERHUB_FILE,
    SERVICE_CLASS_FILE,
    DOCKERHUB_ORG,
    DOCKERHUB_LATEST_TAG,
)


class ModelRegisterer(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def register_from_dockerhub(self):
        data = {"docker_hub": True}
        self.logger.debug(
            "Registering model {0} in the file system".format(self.model_id)
        )
        path = os.path.join(EOS, "dest", self.model_id)
        self.logger.debug(path)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        self.logger.debug(file_name)
        with open(file_name, "w") as f:
            json.dump(data, f)
        current_time = datetime.datetime.now()
        folder_name = current_time.strftime("%Y%m%d%H%M%S")
        path = os.path.join(EOS, "repository", self.model_id)
        if os.path.exists(path):
            shutil.rmtree(path)
        path = os.path.join(path, folder_name)
        os.makedirs(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        file_name = os.path.join(path, SERVICE_CLASS_FILE)
        self.logger.debug("Writing service class pulled_docker {0}".format(file_name))
        with open(file_name, "w") as f:
            f.write("pulled_docker")

    def register_not_from_dockerhub(self):
        data = {"docker_hub": False}
        path = self._model_path(self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        path = self._get_bundle_location(model_id=self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_DOCKERHUB_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)

    def register(self, is_from_dockerhub):
        if is_from_dockerhub:
            self.register_from_dockerhub()
        else:
            self.register_not_from_dockerhub()


class ModelDockerHubFetcher(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.simple_docker = SimpleDocker()

    def is_docker_installed(self):
        return DockerRequirement().is_installed()

    def is_available(self, model_id):
        mp = ModelPuller(model_id=model_id, config_json=self.config_json)
        if mp.is_available_locally():
            return True
        if mp.is_available_in_dockerhub():
            return True
        return False

    def write_apis(self, model_id):
        self.logger.debug("Writing APIs")
        di = PulledDockerImageService(
            model_id=model_id, config_json=self.config_json, preferred_port=None
        )
        di.serve()
        di.close()

    def copy_information(self, model_id):
        fr_file = "/root/eos/dest/{0}/information.json".format(model_id)
        to_file = "{0}/dest/{1}/information.json".format(EOS, model_id)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def copy_metadata(self, model_id):
        fr_file = "/root/eos/dest/{0}/api_schema.json".format(model_id)
        to_file = "{0}/dest/{1}/api_schema.json".format(EOS, model_id)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def copy_status(self, model_id):
        fr_file = "/root/eos/dest/{0}/{1}".format(model_id, STATUS_FILE)
        to_file = "{0}/dest/{1}/{2}".format(EOS, model_id, STATUS_FILE)
        self.simple_docker.cp_from_image(
            img_path=fr_file,
            local_path=to_file,
            org=DOCKERHUB_ORG,
            img=model_id,
            tag=DOCKERHUB_LATEST_TAG,
        )

    def fetch(self, model_id):
        mp = ModelPuller(model_id=model_id, config_json=self.config_json)
        mp.pull()
        mr = ModelRegisterer(model_id=model_id, config_json=self.config_json)
        mr.register(is_from_dockerhub=True)
        self.write_apis(model_id)
        self.copy_information(model_id)
        self.copy_metadata(model_id)
        self.copy_status(model_id)


class ModelFetcher(ErsiliaBase):
    def __init__(
        self,
        config_json=None,
        credentials_json=None,
        overwrite=True,
        repo_path=None,
        mode=None,
        pip=False,
        dockerize=False,
        force_from_github=False,
        force_from_s3=False,
        force_from_dockerhub=False,
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
        self.progress = {}
        self.model_dockerhub_fetcher = ModelDockerHubFetcher(
            config_json=self.config_json
        )
        self.is_docker_installed = self.model_dockerhub_fetcher.is_docker_installed()
        self.force_from_github = force_from_github
        self.force_from_s3 = force_from_s3
        self.force_from_dockerhub = force_from_dockerhub

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
        progress_bar = tqdm(total=8, position=0, leave=True, colour="BLUE")
        start = timer()
        self.model_id = model_id
        self.progress["step0_seconds"] = time.time()
        self._setup_check()
        self.progress["step1_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Checking setup: {0:.3f}s".format(
                self.progress["step1_seconds"] - self.progress["step0_seconds"]
            )
        )
        self._prepare()
        self.progress["step2_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Preparing model: {}s".format(
                self.progress["step2_seconds"] - self.progress["step1_seconds"]
            )
        )
        self._get()
        self.progress["step3_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Getting model: {}s".format(
                self.progress["step3_seconds"] - self.progress["step2_seconds"]
            )
        )
        self._pack()
        self.progress["step4_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Packing model: {}s".format(
                self.progress["step4_seconds"] - self.progress["step3_seconds"]
            )
        )
        self._toolize()
        self.progress["step5_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Checking if model needs to be integrated to a tool: {}s".format(
                self.progress["step5_seconds"] - self.progress["step4_seconds"]
            )
        )
        self._content()
        self.progress["step6_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Getting model card: {}s".format(
                self.progress["step6_seconds"] - self.progress["step5_seconds"]
            )
        )
        self._check()
        self.progress["step7_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Checking that autoservice works: {}s".format(
                self.progress["step7_seconds"] - self.progress["step6_seconds"]
            )
        )
        self._sniff()
        self.progress["step8_seconds"] = time.time()
        progress_bar.update(1)
        tqdm.write(
            "Sniffing model: {}s".format(
                self.progress["step8_seconds"] - self.progress["step7_seconds"]
            )
        )
        self._inform()
        self._success()
        progress_bar.close()
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

    def _decide_if_use_dockerhub(self, model_id):
        if self.repo_path is not None:
            return False
        if self.force_from_dockerhub:
            return True
        if self.force_from_s3:
            return False
        if self.force_from_github:
            return False
        if not self.is_docker_installed:
            self.logger.debug("Docker is not installed in your local")
            return False
        if not self.model_dockerhub_fetcher.is_available(model_id=model_id):
            self.logger.debug("Docker image of this model doesn't seem to be available")
            return False
        return True

    def fetch(self, model_id):
        do_dockerhub = self._decide_if_use_dockerhub(model_id=model_id)
        if do_dockerhub:
            self._fetch_from_dockerhub(model_id=model_id)
        else:
            self._fetch_not_from_dockerhub(model_id=model_id)
