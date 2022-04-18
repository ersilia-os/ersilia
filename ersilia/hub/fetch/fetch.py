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

from ...utils.terminal import run_command
from ...utils.config import Secrets

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

    def _fetchtime(self):
        ts = time.time()
        with open("fetched_models.txt","a") as file:
            file.write(self.model_id)
            file.write(',')
            file.write(str(ts))
            file.write('\n')

    def fetch(self, model_id):
        self.model_id = model_id

        self.progress["step0_seconds"] = time.time()
        self._prepare()

        self.progress["step1_seconds"] = time.time()
        print("step 1 done: {}s".format(self.progress["step1_seconds"] - self.progress["step0_seconds"]))

        self._get()
        self.progress["step2_seconds"] = time.time()
        print("step 2 done: {}s".format(self.progress["step2_seconds"] - self.progress["step1_seconds"]))

        self._pack()
        self.progress["step3_seconds"] = time.time()
        print("step 3 done: {}s".format(self.progress["step3_seconds"]-self.progress["step2_seconds"]))

        self._toolize()
        self.progress["step4_seconds"] = time.time()
        print("step 4 done: {}s".format(self.progress["step4_seconds"] - self.progress["step3_seconds"]))
        
        self._content()
        self.progress["step5_seconds"] = time.time()
        print("step 5 done: {}s".format(self.progress["step5_seconds"] - self.progress["step4_seconds"]))

        self._check()
        self.progress["step6_seconds"] = time.time()
        print("step 6 done: {}s".format(self.progress["step6_seconds"] - self.progress["step5_seconds"]))
        
        self._sniff()
        self.progress["step7_seconds"] = time.time()
        print("step 7 done: {}s".format(self.progress["step7_seconds"] - self.progress["step6_seconds"]))

        self._success()
        print("Fetching {0} done in time: {1}s".format(model_id, abs(self.progress["step7_seconds"]-self.progress["step0_seconds"])))
        for i in reversed(range(1, len(self.progress))):
            self.progress["step{0}_seconds".format(i)] -= self.progress["step{0}_seconds".format(i-1)]

        with open("METADATA.json", "w") as outfile:
            json.dump({model_id: {"progress": self.progress}}, outfile)

        content = base64.b64encode(json.dumps(self.progress).encode('utf-8'))
        data = json.dumps({"message":"[update] METADATA from Ersilia","content": str(content)})
        token = "GITHUB_ACCESS_TOKEN"
        cmd = 'curl -X PUT -H "Accept: application/vnd.github.v3+json" -H "Authorization: token {0}" https://api.github.com/repos/ersilia-os/{1}/contents/METADATA.json -d {2}'.format(token, self.model_id, data)
        run_command(cmd)

        # print("Progress times for each step in seconds", self.progress)
