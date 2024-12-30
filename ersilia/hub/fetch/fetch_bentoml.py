"""Fetch model from the Ersilia Model Hub using BentoML."""

import json
import os
from timeit import default_timer as timer
from datetime import timedelta

from .actions.template_resolver import TemplateResolver
from .actions.setup import SetupChecker
from .actions.prepare import ModelPreparer
from .actions.get import ModelGetter
from .actions.lake import LakeGetter
from .actions.pack_bentoml import ModelPacker
from .actions.toolize import ModelToolizer
from .actions.content import CardGetter
from .actions.check import ModelChecker
from .actions.sniff_bentoml import ModelSniffer
from .actions.inform import ModelInformer
from .register.register import ModelRegisterer

from ... import ErsiliaBase

from . import STATUS_FILE, DONE_TAG


class ModelFetcherFromBentoML(ErsiliaBase):
    """
    ModelFetcherFromBentoML is responsible for fetching models from the Ersilia Model Hub using BentoML.

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

    Examples
    --------
    .. code-block:: python

        fetcher = ModelFetcherFromBentoML(config_json=config)
        fetcher.fetch(model_id="eosxxxx")
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
            force_from_github=self.force_from_github,
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

    def _fetch_not_from_dockerhub(self, model_id: str):
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

    def _fetch(self, model_id: str):
        self.logger.debug("Starting fetching procedure")
        if self.overwrite is None:
            self.overwrite = True
        self._fetch_not_from_dockerhub(model_id=model_id)

    def seems_installable(self, model_id: str) -> bool:
        """
        Check if the model is installable with BentoML.

        Parameters
        ----------
        model_id : str
            The ID of the model to be checked.

        Returns
        -------
        bool
            True if the model is installable with BentoML, False otherwise.
        """
        tr = TemplateResolver(
            model_id=model_id, repo_path=self.repo_path, config_json=self.config_json
        )
        self.logger.debug("Checking if the model is installable with BentoML")
        return tr.is_bentoml()

    def fetch(self, model_id: str):
        """
        Fetch the model from the Ersilia Model Hub.

        This method initiates the fetching process for the model.

        Parameters
        ----------
        model_id : str
            The ID of the model to be fetched.

        Examples
        --------
        .. code-block:: python

            fetcher = ModelFetcherFromBentoML(config_json=config)
            fetcher.fetch(model_id="eosxxxx")
        """
        self._fetch(model_id)
