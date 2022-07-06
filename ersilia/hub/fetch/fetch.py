"""Fetch model from the Ersilia Model Hub"""

import json
import os
import time

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

    def _fetchtime(self):
        ts = time.time()
        with open("fetched_models.txt", "a") as file:
            file.write(self.model_id)
            file.write(",")
            file.write(str(ts))
            file.write("\n")

    def fetch(self, model_id):
        self.model_id = model_id
        self._prepare()
        self._get()
        self._pack()
        self._toolize()
        self._content()
        self._check()
        self._sniff()
        self._success()
        # self._fetchtime()
        logger.info("Fetching {0} done successfully".format(model_id))
