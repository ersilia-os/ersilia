"""Fetch model from the Ersilia Model Hub"""

import json
import os
import time

from numpy import subtract

from ... import ErsiliaBase
from ... import logger
from .actions.prepare import ModelPreparer
from .actions.get import ModelGetter
from .actions.lake import LakeGetter
from .actions.pack import ModelPacker
from .actions.toolize import ModelToolizer
from .actions.content import CardGetter
from .actions.check import ModelChecker
from .actions.sniff import ModelSniffer

from . import STATUS_FILE, DONE_TAG


class ModelFetcher(ErsiliaBase):
    def __init__(
        self,
        config_json=None,
        credentials_json=None,
        overwrite=True,
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
        if self.mode == "docker":
            self.logger.debug("When packing mode is docker, dockerization is mandatory")
            dockerize = True
        self.do_docker = dockerize
        self.progress = {}

    def _prepare(self):
        mp = ModelPreparer(
            model_id=self.model_id,
            overwrite=self.overwrite,
            config_json=self.config_json,
        )
        mp.prepare()

    def _get(self):
        mg = ModelGetter(self.model_id, self.config_json)
        mg.get()

    def _lake(self):
        ml = LakeGetter(self.model_id, self.config_json)
        ml.get()

    def _pack(self):
        mp = ModelPacker(self.model_id, self.mode, self.config_json)
        mp.pack()

    def _toolize(self):
        mt = ModelToolizer(self.model_id, self.config_json)
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

    def _success(self):
        done = {DONE_TAG: True}
        status_file = os.path.join(self._dest_dir, self.model_id, STATUS_FILE)
        with open(status_file, "w") as f:
            json.dump(done, f, indent=4)

    def fetch(self, model_id):
        self.model_id = model_id

        self.progress["step0_seconds"] = time.time()
        self._prepare()
        # step1_seconds = time.time() - start_seconds
        self.progress["step1_seconds"] = time.time()
        print("step 1 done")

        # step2_seconds = time.time() - step1_seconds
        self._get()
        print("step 2 done")
        self.progress["step2_seconds"] = time.time()

        # step3_seconds = time.time() - step2_seconds
        self._pack()
        print("step 3 done")
        self.progress["step3_seconds"] = time.time()

        # step4_seconds = time.time() - step3_seconds
        self._toolize()
        print("step 4 done")
        self.progress["step4_seconds"] = time.time()

        # step5_seconds = time.time() - step4_seconds
        self._content()
        print("step 5 done")
        self.progress["step5_seconds"] = time.time()

        # step6_seconds = time.time() - step5_seconds
        self._check()
        print("step 6 done")
        self.progress["step6_seconds"] = time.time()

        # step7_seconds = time.time() - step6_seconds
        self._sniff()
        print("step 7 done")
        # final_seconds = time.time() - step7_seconds
        self.progress["step7_seconds"] = time.time()

        self._success()
        print("Fetching {0} done in time: {1}s".format(model_id, abs(self.progress["step7_seconds"]-self.progress["step0_seconds"])))
        for i in reversed(range(1, len(self.progress))):
            self.progress["step{0}_seconds".format(i)] -= self.progress["step{0}_seconds".format(i-1)]

        with open("{0}_progress.json".format(model_id), "w") as outfile:
            json.dump(self.progress, outfile)
        print("Progress times for each step in seconds", self.progress)
