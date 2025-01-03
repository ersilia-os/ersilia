"""Fetch model from the Ersilia Model Hub using FastAPI."""

import json
import os
from datetime import timedelta
from timeit import default_timer as timer

from ... import ErsiliaBase
from . import DONE_TAG, STATUS_FILE
from .actions.check import ModelChecker
from .actions.content import CardGetter
from .actions.get import ModelGetter
from .actions.inform import ModelInformer
from .actions.pack_fastapi import ModelPacker
from .actions.prepare import ModelPreparer
from .actions.setup import SetupChecker
from .actions.sniff_fastapi import ModelSniffer
from .actions.template_resolver import TemplateResolver
from .register.register import ModelRegisterer


class ModelFetcherFromFastAPI(ErsiliaBase):
    """
    ModelFetcherFromFastAPI is responsible for fetching models from the Ersilia Model Hub using FastAPI.

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
    force_from_github : bool, optional
        Whether to force fetching from GitHub.
    force_from_s3 : bool, optional
        Whether to force fetching from S3.

    Examples
    --------
    .. code-block:: python

        fetcher = ModelFetcherFromFastAPI(
            config_json=config
        )
        fetcher.fetch(model_id="eosxxxx")
    """

    def __init__(
        self,
        config_json: dict = None,
        credentials_json: dict = None,
        overwrite: bool = None,
        repo_path: str = None,
        mode: str = None,
        force_from_github: bool = False,
        force_from_s3: bool = False,
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
            force_from_github=self.force_from_github,
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

    def _fetch(self, model_id: str):
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

    def seems_installable(self, model_id: str) -> bool:
        """
        Check if the model is installable with FastAPI.

        Parameters
        ----------
        model_id : str
            The ID of the model to be checked.

        Returns
        -------
        bool
            True if the model is installable with FastAPI, False otherwise.
        """
        tr = TemplateResolver(
            model_id=model_id, repo_path=self.repo_path, config_json=self.config_json
        )
        return tr.is_fastapi()

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

            fetcher = ModelFetcherFromFastAPI(
                config_json=config
            )
            fetcher.fetch(model_id="eosxxxx")
        """
        self.logger.debug("Fetching from FastAPI...")
        self._fetch(model_id=model_id)
