"""Fetch model from the Ersilia Model Hub using FastAPI."""

import os
import json

from timeit import default_timer as timer
from datetime import timedelta

from .actions.template_resolver import TemplateResolver
from .actions.setup import SetupChecker
from .actions.prepare import ModelPreparer
from .actions.get import ModelGetter
from .actions.pack_fastapi import ModelPacker
from .actions.content import CardGetter
from .actions.check import ModelChecker
from .actions.sniff_fastapi import ModelSniffer
from .actions.inform import ModelInformer
from .register.register import ModelRegisterer

from ... import ErsiliaBase

from . import STATUS_FILE, DONE_TAG


class ModelFetcherFromFastAPI(ErsiliaBase):
    def __init__(
        self,
        config_json=None,
        credentials_json=None,
        overwrite=None,
        repo_path=None,
        mode=None,
        force_from_github=False,
        force_from_s3=False,
    ):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.repo_path = repo_path
        self.overwrite = overwrite
        self.mode = mode
        self.force_from_github = force_from_github
        self.force_from_s3 = force_from_s3

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

    def _pack(self):
        mp = ModelPacker(
            model_id=self.model_id, mode=self.mode, config_json=self.config_json
        )
        mp.pack()

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

    def _fetch(self, model_id):
        start = timer()
        self.model_id = model_id
        self._setup_check()
        self._prepare()
        self._get()
        self._pack()
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

    def seems_installable(self, model_id):
        tr = TemplateResolver(
            model_id=model_id, repo_path=self.repo_path, config_json=self.config_json
        )
        return tr.is_fastapi()

    def fetch(self, model_id):
        self.logger.debug("Fetching from FastAPI...")
        self._fetch(model_id=model_id)
