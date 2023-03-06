"""Fetch model from the Ersilia Model Hub"""

import json
import os
from timeit import default_timer as timer
from datetime import timedelta
import time
from tqdm import tqdm

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

from . import STATUS_FILE, DONE_TAG
from ...default import INFORMATION_FILE


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

    def fetch(self, model_id):
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
        print("Fetching {0} done in time: {1}s".format(model_id, abs(elapsed_time)))

        self.logger.info(
            "Fetching {0} done successfully: {1}".format(model_id, elapsed_time)
        )
