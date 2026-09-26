import json
import os
import re

from .. import ErsiliaBase, throw_ersilia_exception
from ..default import DOCKER_INFO_FILE
from ..hub.content.slug import Slug
from ..hub.fetch import DONE_TAG, STATUS_FILE
from ..utils.exceptions_utils.exceptions import InvalidModelIdentifierError
from ..utils.paths import get_metadata_from_base_dir


def _suggest_model(text):
    # The closest models in the Hub (up to three), e.g. "eos3b5e
    # (molecular-weight)", or None (also when the Hub cannot be reached).
    import difflib

    try:
        from ..db.hubdata.interfaces import JsonModelsInterface

        models = JsonModelsInterface().items_all()
    except Exception:
        return None
    names = {}
    for m in models:
        identifier, slug = m.get("Identifier"), m.get("Slug")
        if identifier:
            label = f"{identifier} ({slug})" if slug else identifier
            names[identifier] = label
            if slug:
                names[slug] = label
    matches = difflib.get_close_matches(text.lower(), list(names), n=6, cutoff=0.7)
    labels = list(dict.fromkeys(names[m] for m in matches))[:3]
    if not labels:
        return None
    if len(labels) == 1:
        return labels[0]
    return ", ".join(labels[:-1]) + " or " + labels[-1]


class ModelBase(ErsiliaBase):
    """
    Base class for managing models.

    This class provides foundational functionality for handling models, including initialization,
    validation, and checking local availability.

    Parameters
    ----------
    model_id_or_slug : str, optional
        The model identifier or slug, by default None.
    repo_path : str, optional
        The repository path, by default None.
    config_json : dict, optional
        Configuration in JSON format, by default None.
    """

    @throw_ersilia_exception()
    def __init__(self, model_id_or_slug=None, repo_path=None, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        if model_id_or_slug is None and repo_path is None:
            raise Exception
        if model_id_or_slug is not None and repo_path is not None:
            raise Exception
        if model_id_or_slug is not None:
            model_id_or_slug = model_id_or_slug.strip()
            # Identifiers are lower case (e.g. EOS3B5E is eos3b5e).
            if re.fullmatch(r"(?i)eos[0-9][a-z0-9]{3}", model_id_or_slug):
                model_id_or_slug = model_id_or_slug.lower()
            self.text = model_id_or_slug
            slugger = Slug()
            if slugger.is_slug(model_id_or_slug):
                self.slug = model_id_or_slug
                self.model_id = slugger.encode(self.slug)
            else:
                self.model_id = model_id_or_slug
                self.slug = slugger.decode(self.model_id)
            if not self.is_valid():
                raise InvalidModelIdentifierError(
                    model=self.text, suggestion=_suggest_model(self.text)
                )

        if repo_path is not None:
            self.logger.debug(f"Repo path specified: {repo_path}")
            expanded_path = os.path.expanduser(repo_path)
            abspath = os.path.abspath(expanded_path)
            self.logger.debug(f"Absolute path: {abspath}")
            # Check if path actually exists
            if not os.path.exists(abspath):
                raise FileNotFoundError(
                    "Model directory does not exist at the provided path. Please check the path and try again."
                )
            self.logger.debug(f"Path exists: {abspath}")
            self.text = self._get_model_id_from_path(repo_path)
            self.logger.debug(f"Model ID from path: {self.text}")
            self.model_id = self.text
            slug = self._get_slug_if_available(repo_path)
            if slug is None:
                self.slug = "my-model"
            else:
                self.slug = slug
            self.logger.debug(f"Slug from path: {self.slug}")

    def _get_model_id_from_path(self, repo_path):
        return os.path.basename(os.path.abspath(repo_path)).rstrip("/")

    def _get_slug_if_available(self, repo_path):
        try:
            data = get_metadata_from_base_dir(repo_path)
        except FileNotFoundError:
            return None
        slug = data["Slug"]
        if slug == "":
            return None
        else:
            return slug

    def is_valid(self):
        """
        Check if the model identifier and slug are valid.

        Returns
        -------
        bool
            True if the model identifier and slug are valid, False otherwise.
        """
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
            self._dest_dir, self.model_id, DOCKER_INFO_FILE
        )
        if not os.path.exists(from_dockerhub_file):
            return False
        else:
            return True

    def is_available_locally(self):
        """
        Check if the model is available locally either from the status file
        or from DockerHub.

        Returns
        -------
        bool
            True if the model is available locally, False otherwise.
        """
        bs = self._is_available_locally_from_status()
        bd = self._is_available_locally_from_dockerhub()
        if bs or bd:
            return True
        else:
            return False

    def was_fetched_from_dockerhub(self):
        """
        Check if the model was fetched from DockerHub by reading the DockerHub file.

        Returns
        -------
        bool
            True if the model was fetched from DockerHub, False otherwise.
        """
        from_dockerhub_file = os.path.join(
            self._dest_dir, self.model_id, DOCKER_INFO_FILE
        )
        if not os.path.exists(from_dockerhub_file):
            return False
        with open(from_dockerhub_file, "r") as f:
            data = json.load(f)
            return data["docker_hub"]
