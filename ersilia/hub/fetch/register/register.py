import datetime
import json
import os
import shutil

import validators

from .... import EOS, ErsiliaBase, throw_ersilia_exception
from ....default import (
    DOCKER_INFO_FILE,
    IS_FETCHED_FROM_HOSTED_FILE,
    SERVICE_CLASS_FILE,
)
from ....utils.exceptions_utils.fetch_exceptions import InvalidUrlError
from ...fetch import ModelURLResolver


class ModelRegisterer(ErsiliaBase):
    """
    ModelRegisterer is responsible for registering models from various sources.

    Parameters
    ----------
    model_id : str
        The ID of the model to be registered.
    config_json : dict
        Configuration settings for the registerer.

    Examples
    --------
    .. code-block:: python

        registerer = ModelRegisterer(
            model_id="eosxxxx", config_json=config
        )
        await registerer.register(is_from_dockerhub=True)
    """

    def __init__(self, model_id: str, config_json: dict):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def register_from_dockerhub(self, **kwargs):
        """
        Register the model from DockerHub.

        This method registers the model in the file system indicating it was fetched from DockerHub.
        """
        if "img_tag" in kwargs:
            img_tag = kwargs["img_tag"]
            data = {"docker_hub": True, "tag": img_tag}
        else:
            data = {"docker_hub": True}
        self.logger.debug(
            "Registering model {0} in the file system".format(self.model_id)
        )
        path = os.path.join(EOS, "dest", self.model_id)
        self.logger.debug(path)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        file_name = os.path.join(path, DOCKER_INFO_FILE)
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
        file_name = os.path.join(path, DOCKER_INFO_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        file_name = os.path.join(path, SERVICE_CLASS_FILE)
        self.logger.debug("Writing service class pulled_docker {0}".format(file_name))
        with open(file_name, "w") as f:
            f.write("pulled_docker")

    def register_not_from_dockerhub(self):
        """
        Register the model indicating it was not fetched from DockerHub.

        This method registers the model in the file system indicating it was not fetched from DockerHub.
        """
        data = {"docker_hub": False}
        path = self._model_path(self.model_id)
        file_name = os.path.join(path, DOCKER_INFO_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        path = self._get_bundle_location(model_id=self.model_id)
        file_name = os.path.join(path, DOCKER_INFO_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)

    def _resolve_url(self):
        """
        Resolve the URL for the hosted model.

        This method resolves a valid URL for the hosted model.

        Returns
        -------
        str or None
            The resolved URL if valid, otherwise None.
        """
        mdl_url_resolver = ModelURLResolver(
            model_id=self.model_id, config_json=self.config_json
        )
        is_valid_url, url = mdl_url_resolver.resolve_valid_hosted_model_url(
            self.model_id
        )
        if is_valid_url:
            return url
        else:
            return None

    @throw_ersilia_exception()
    def register_from_hosted(self, url: str = None):
        """
        Register the model from a hosted URL.

        This method registers the model in the file system indicating it was fetched from a hosted URL.

        Parameters
        ----------
        url : str, optional
            The URL from which the model is hosted. If not provided, it will be resolved automatically.
        """
        if url is None:
            url = self._resolve_url()
        else:
            if not validators.url(url):
                raise InvalidUrlError(url)
        data = {"hosted": True, "url": url}
        self.logger.debug(
            "Registering model {0} in the file system".format(self.model_id)
        )
        path = os.path.join(EOS, "dest", self.model_id)
        self.logger.debug(path)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
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
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        file_name = os.path.join(path, SERVICE_CLASS_FILE)
        self.logger.debug("Writing service class hosted {0}".format(file_name))
        with open(file_name, "w") as f:
            f.write("hosted")

    def register_not_from_hosted(self):
        """
        Register the model indicating it was not fetched from a hosted URL.

        This method registers the model in the file system indicating it was not fetched from a hosted URL.
        """
        data = {"hosted": False}
        path = self._model_path(self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)
        path = self._get_bundle_location(model_id=self.model_id)
        file_name = os.path.join(path, IS_FETCHED_FROM_HOSTED_FILE)
        with open(file_name, "w") as f:
            json.dump(data, f)

    async def register(
        self, is_from_dockerhub: bool = False, is_from_hosted: bool = False, **kwargs
    ):
        """
        Register the model based on its source.

        This method registers the model in the file system based on whether it was fetched from DockerHub or a hosted URL.

        Parameters
        ----------
        is_from_dockerhub : bool, optional
            Indicates if the model is from DockerHub.
        is_from_hosted : bool, optional
            Indicates if the model is from a hosted URL.

        Raises
        ------
        ValueError
            If both is_from_dockerhub and is_from_hosted are True.

        Examples
        --------
        .. code-block:: python

            registerer = ModelRegisterer(
                model_id="eosxxxx", config_json=config
            )
            await registerer.register(is_from_dockerhub=True)
        """
        if is_from_dockerhub and is_from_hosted:
            raise ValueError("Model cannot be from both DockerHub and hosted")
        elif is_from_dockerhub and not is_from_hosted:
            self.register_from_dockerhub(**kwargs)
            self.register_not_from_hosted()
        elif not is_from_dockerhub and is_from_hosted:
            self.register_from_hosted()
            self.register_not_from_dockerhub()
        else:
            self.register_not_from_dockerhub()
            self.register_not_from_hosted()
