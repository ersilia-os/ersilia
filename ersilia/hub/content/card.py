import json
import os

import requests
import yaml

from ... import ErsiliaBase
from ...auth.auth import Auth
from ...db.hubdata.interfaces import JsonModelsInterface
from ...default import (
    CARD_FILE,
    INFORMATION_FILE,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
)
from ...utils.logging import make_temp_dir
from ...utils.paths import get_metadata_from_base_dir
from ...utils.terminal import run_command
from .base_information import BaseInformation


class RepoMetadataFile(ErsiliaBase):
    """
    Class to handle the metadata file of a model repository.

    This class provides methods to get the URL of the metadata file from GitHub,
    read information from the metadata file, and write information to the metadata file.

    Parameters
    ----------
    model_id : str, optional
        The model identifier.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, model_id=None, config_json=None):
        self.model_id = model_id
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def _github_json_url(self, org=None, branch=None):
        if org is None:
            org = "ersilia-os"
        if branch is None:
            branch = "main"
        return "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
            org, self.model_id, branch, METADATA_JSON_FILE
        )

    def _github_yaml_url(self, org=None, branch=None):
        if org is None:
            org = "ersilia-os"
        if branch is None:
            branch = "main"
        return "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
            org, self.model_id, branch, METADATA_YAML_FILE
        )

    def _get_file_content_from_github(self, org, branch):
        json_url = self._github_json_url(org, branch)
        r = requests.get(json_url)
        if r.status_code == 404:
            yaml_url = self._github_yaml_url(org, branch)
            r = requests.get(yaml_url)
            if r.status_code == 404:
                return None
            else:
                return yaml.safe_load(r.content)
        else:
            return json.loads(r.content)

    def get_json_or_yaml_file(self, org: str = None, branch: str = None) -> dict:
        """
        Get the metadata file from GitHub in JSON or YAML format. JSON format typically used for bentoml
        packed models and YAML format typically used for ersilia pack models.

        Parameters
        ----------
        org : str, optional
            The GitHub organization.
        branch : str, optional
            The GitHub branch.

        Returns
        -------
        dict
            The metadata file content.
        """
        return self._get_file_content_from_github(org, branch)

    def read_information(
        self,
        org: str = None,
        branch: str = None,
        json_or_yaml_path: str = None,
    ) -> BaseInformation:
        """
        Read information from the metadata file.

        Parameters
        ----------
        org : str, optional
            The GitHub organization.
        branch : str, optional
            The GitHub branch.
        json_or_yaml_path : str, optional
            The path to the JSON or YAML file.

        Returns
        -------
        BaseInformation
            The base information read from the metadata file.
        """
        if json_or_yaml_path is None:
            data = self.get_json_or_yaml_file(org=org, branch=branch)
        else:
            with open(json_or_yaml_path, "r") as f:
                if json_or_yaml_path.endswith(".json"):
                    data = json.load(f)
                elif json_or_yaml_path.endswith(".yaml") or json_or_yaml_path.endswith(
                    ".yml"
                ):
                    data = yaml.safe_load(f)
                else:
                    raise ValueError("File format not supported")
        bi = BaseInformation(config_json=self.config_json)
        bi.from_dict(data)
        return bi

    def write_information(
        self, data: BaseInformation, json_or_yaml_path: str = None
    ) -> dict:
        """
        Write information to the metadata file.

        Parameters
        ----------
        data : BaseInformation
            The base information to write.
        json_or_yaml_path : str, optional
            The path to the JSON or YAML file.

        Returns
        -------
        dict
            The written data.
        """
        data = data.as_dict()
        if json_or_yaml_path.endswith(".json"):
            path = json_or_yaml_path
            with open(path, "w") as f:
                json.dump(data, f, indent=4)
        elif json_or_yaml_path.endswith(".yml") or json_or_yaml_path.endswith(".yaml"):
            path = json_or_yaml_path
            with open(path, "w") as f:
                yaml.dump(data, f)
        else:
            raise ValueError("File format not supported")
        return data


class MetadataCard(ErsiliaBase):
    """
    Class to handle the metadata card of a model.

    This class provides methods to get the metadata card of a model from the local repository.

    Parameters
    ----------
    config_json : dict
        Configuration settings in JSON format.
    """

    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.logger.debug("Initializing MetadataCard")

    def get(self, model_id: str = None) -> dict:
        """
        Get the metadata card of a model from the local repository.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.

        Returns
        -------
        dict
            The metadata card of the model.
        """
        self.logger.debug("Getting metadata card for model_id: {0}".format(model_id))
        if model_id is not None:
            dest_dir = self._model_path(model_id=model_id)
            self.logger.debug("Trying to get metadata from: {0}".format(dest_dir))
            try:
                data = get_metadata_from_base_dir(dest_dir)
            except FileNotFoundError:
                return None
            return data
        else:
            return


class ReadmeCard(ErsiliaBase):
    """
    Class to handle the README card of a model.

    This class provides methods to get the README card of a model from GitHub.

    Parameters
    ----------
    config_json : dict
        Configuration settings in JSON format.
    """

    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _raw_readme_url(self, model_id):
        url = (
            "https://raw.githubusercontent.com/ersilia-os/{0}/master/README.md".format(
                model_id
            )
        )
        return url

    def _gh_view(self, model_id):
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "view.md")
        cmd = "gh repo view {0}/{1} > {2}".format("ersilia-os", model_id, tmp_file)
        run_command(cmd)
        with open(tmp_file, "r") as f:
            text = f.read()
        return text

    def _title(self, lines):
        """Title is determined by the first main header in markdown"""
        for l in lines:
            if l[:2] == "# ":
                s = l[2:].strip()
                return s

    def _description(self, lines):
        """Description is what comes after the title and before the next header"""
        text = "\n".join(lines)
        return text.split("# ")[1].split("\n")[1].split("#")[0].strip()

    def _model_github_url(self, model_id):
        return "https://github.com/ersilia-os/{0}".format(model_id)

    def parse(self, model_id):
        """
        Parse the model information from the README file.

        Parameters
        ----------
        model_id : str
            The model identifier.
        """
        readme = os.path.join(self._dest_dir, model_id, "README.md")
        if os.path.exists(readme):
            with open(readme, "r") as f:
                text = f.read()
        else:
            if Auth().is_contributor():
                text = self._gh_view(model_id)
                if not text:
                    return None
                text = "--".join(text.split("--")[1:])
            else:
                r = requests.get(self._raw_readme_url(model_id))
                if r.status_code != 200:
                    return None
                text = r.text
        lines = text.split(os.linesep)
        results = {
            "model_id": model_id,
            "title": self._title(lines),
            "description": self._description(lines),
            "github_url": self._model_github_url(model_id),
        }
        return results

    def get(self, model_id: str = None, slug: str = None) -> dict:
        """
        Get the README card of a model from GitHub.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.
        slug : str, optional
            The slug of the model.

        Returns
        -------
        dict
            The README card of the model.
        """
        if model_id:
            return self.parse(model_id)
        else:
            return None


class LocalCard(ErsiliaBase):
    """
    Class to handle the local card of a model.

    This class provides methods to get the local card of a model from the local repository stored in eos directory.

    Parameters
    ----------
    config_json : dict
        Configuration settings in JSON format.
    """

    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _load_data(self, model_id):
        """
        Loads the JSON data from the model's information file.
        """
        model_path = self._model_path(model_id)
        info_file = os.path.join(model_path, INFORMATION_FILE)
        if os.path.exists(info_file):
            card_path = info_file
        else:
            card_path = os.path.join(model_path, CARD_FILE)
        if os.path.exists(card_path):
            with open(card_path, "r") as f:
                card = json.load(f)
            return card
        else:
            return None

    def get(self, model_id: str = None, slug: str = None) -> dict:
        """
        Get the local card of a model from the local repository.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.
        slug : str, optional
            The slug of the model.

        Returns
        -------
        dict
            The local card of the model.
        """
        if model_id:
            card = self._load_data(model_id)
            if card is None:
                mc = MetadataCard(config_json=self.config_json)
                card = mc.get(model_id)
            return card
        else:
            return


class S3JsonCard(JsonModelsInterface):
    """
    Class to handle the S3 JSON card of a model.

    This class provides methods to get the S3 JSON card of a model from the S3 bucket.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, config_json=None):
        JsonModelsInterface.__init__(self, config_json=config_json)

    def get_card_by_model_id(self, model_id):
        """
        Get the card information by model identifier.

        Parameters
        ----------
        model_id : str
            The model identifier.

        Returns
        -------
        dict
            The card information.
        """
        all_models = self.items_all()
        for model in all_models:
            if model["Identifier"] == model_id:
                return model

    def get_card_by_slug(self, slug):
        """
        Get the card information by model slug.

        Parameters
        ----------
        slug : str
            The model slug.

        Returns
        -------
        dict
            The card information.
        """
        all_models = self.items_all()
        for model in all_models:
            if model["Slug"] == slug:
                return model

    def get(self, model_id: str = None, slug: str = None) -> dict:
        """
        Get the S3 JSON card of a model from the S3 bucket by model id or slug.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.
        slug : str, optional
            The slug of the model.

        Returns
        -------
        dict
            The S3 JSON card of the model.
        """
        if model_id is not None:
            return self.get_card_by_model_id(model_id)
        elif slug is not None:
            return self.get_card_by_slug(slug)
        else:
            raise ValueError("Either model_id or slug must be provided")


class ModelCard(object):
    """
    Class to handle the model card.

    This class provides methods to get the model card from various sources,
    such as local repository, metadata card, S3 JSON card, and README card.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, config_json=None):
        self.config_json = config_json

    def _get(self, model_id, slug):
        lc = LocalCard(config_json=self.config_json)
        card = lc.get(model_id, slug)
        if card is not None:
            return card
        mc = MetadataCard(config_json=self.config_json)
        card = mc.get(model_id)
        if card is not None:
            return card
        jc = S3JsonCard(config_json=self.config_json)
        card = jc.get(model_id, slug)
        if card is not None:
            return card
        rc = ReadmeCard(config_json=self.config_json)
        card = rc.get(model_id, slug)
        if card is not None:
            return card

    def get(
        self,
        model_id: str = None,
        slug: str = None,
        as_json: bool = False,
    ) -> dict:
        """
        Get the model card from various sources.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.
        slug : str, optional
            The slug of the model.
        as_json : bool, optional
            Whether to return the model card in JSON format.

        Returns
        -------
        dict
            The model card.
        """
        card = self._get(model_id, slug)
        if card is None:
            return
        if as_json:
            return json.dumps(card, indent=4)
        else:
            return card
