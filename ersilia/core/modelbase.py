import os
import json

from .. import ErsiliaBase
from ..hub.content.slug import Slug
from ..hub.fetch import STATUS_FILE, DONE_TAG
from ..default import IS_FETCHED_FROM_DOCKERHUB_FILE

from ..utils.exceptions_utils.exceptions import InvalidModelIdentifierError
from .. import throw_ersilia_exception


class ModelBase(ErsiliaBase):
    """Base class of a Model."""

    @throw_ersilia_exception
    def __init__(self, model_id_or_slug=None, repo_path=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        if model_id_or_slug is None and repo_path is None:
            raise Exception
        if model_id_or_slug is not None and repo_path is not None:
            raise Exception
        if model_id_or_slug is not None:
            self.text = model_id_or_slug
            slugger = Slug()

            if slugger.is_slug(model_id_or_slug):
                self.slug = model_id_or_slug
                self.model_id = slugger.encode(self.slug)
            else:
                self.model_id = model_id_or_slug
                self.slug = slugger.decode(self.model_id)
            if not self.is_valid():
                raise InvalidModelIdentifierError(model=self.text)

        if repo_path is not None:
            self.logger.debug("Repo path specified: {0}".format(repo_path))
            self.logger.debug("Absolute path: {0}".format(os.path.abspath(repo_path)))
            self.text = self._get_model_id_from_path(repo_path)
            self.model_id = self.text
            slug = self._get_slug_if_available(repo_path)
            if slug is None:
                self.slug = "my-model"
            else:
                self.slug = slug

    def _get_model_id_from_path(self, repo_path):
        return os.path.basename(os.path.abspath(repo_path)).rstrip("/")

    def _get_slug_if_available(self, repo_path):
        metadata_json = os.path.join(repo_path, "metadata.json")
        if os.path.exists(metadata_json):
            with open(metadata_json, "r") as f:
                data = json.load(f)
            slug = data["Slug"]
            if slug == "":
                return None
            else:
                return slug
        else:
            return None

    def is_valid(self):
        if self.model_id is None or self.slug is None:
            return False
        else:
            return True

    def _is_available_locally_from_status(self):
        fetch_status_file = os.path.join(self._dest_dir, self.model_id, STATUS_FILE)
        if not os.path.exists(fetch_status_file):
            self.logger.debug("No status file exists")
            is_fetched = False
        else:
            with open(fetch_status_file, "r") as f:
                status = json.load(f)
            is_fetched = status[DONE_TAG]
        self.logger.debug("Is fetched: {0}".format(is_fetched))
        return is_fetched

    def _is_available_locally_from_dockerhub(self):
        from_dockerhub_file = os.path.join(
            self._dest_dir, self.model_id, IS_FETCHED_FROM_DOCKERHUB_FILE
        )
        if not os.path.exists(from_dockerhub_file):
            return False
        else:
            return True

    def is_available_locally(self):
        bs = self._is_available_locally_from_status()
        bd = self._is_available_locally_from_dockerhub()
        if bs or bd:
            return True
        else:
            return False

    def was_fetched_from_dockerhub(self):
        from_dockerhub_file = os.path.join(
            self._dest_dir, self.model_id, IS_FETCHED_FROM_DOCKERHUB_FILE
        )
        if not os.path.exists(from_dockerhub_file):
            return False
        with open(from_dockerhub_file, "r") as f:
            data = json.load(f)
            return data["docker_hub"]
