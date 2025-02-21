import json
import os

import requests
import yaml

from ... import ErsiliaBase
from ...db.hubdata.interfaces import JsonModelsInterface
from .base_information import BaseInformation

try:
    from isaura.core.hdf5 import Hdf5Explorer
except:
    Hdf5Explorer = None

from ...default import (
    CARD_FILE,
    INFORMATION_FILE,
    METADATA_JSON_FILE,
    METADATA_YAML_FILE,
)
from ...utils.paths import get_metadata_from_base_dir


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
                elif json_or_yaml_path.endswith(".yaml"):
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
        elif json_or_yaml_path.endswith(".yaml"):
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

    def get(self, model_id: str = None, slug: str = None) -> dict:
        """
        Get the metadata card of a model from the local repository.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.
        slug : str, optional
            The slug of the model.

        Returns
        -------
        dict
            The metadata card of the model.
        """
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


class ReadmeCardVersion0(ErsiliaBase):
    """
    Parser for the old README format.
    """

    def __init__(self, config_json):
        super().__init__(config_json=config_json)

    def _parse_title(self, lines):
        for line in lines:
            if line.startswith("# "):
                return line[2:].strip()
        return "Unknown Title"

    def _parse_description(self, lines):
        text = "\n".join(lines)
        try:
            return text.split("# ")[1].split("\n", 1)[1].split("##")[0].strip()
        except IndexError:
            return "No description found."

    def parse_text(self, text, model_id):
        """
        Parses the text from an old-format README file and returns a dictionary of metadata.
        Parameters:
            text (str): The README file content.
            model_id (str): The model identifier.
        Returns:
            dict: Parsed metadata.
        """
        lines = text.splitlines()
        results = {
            "model_id": model_id,
            "title": self._parse_title(lines),
            "description": self._parse_description(lines),
            "github_url": "https://github.com/ersilia-os/{}".format(model_id),
        }
        return results


class ReadmeCardVersion1(ErsiliaBase):
    """
    Parser for the new README format.
    Extracts sections (e.g., Description, Identifiers, Domain, etc.) based on header markers.
    """

    def __init__(self, config_json):
        super().__init__(config_json=config_json)

    def _extract_section(self, text, header, next_marker="\n## "):
        """
        Extracts text from a section starting with the given header.
        Extraction stops at the next header marker.
        """
        if header in text:
            section = text.split(header, 1)[1]
            if next_marker in section:
                section = section.split(next_marker, 1)[0]
            return section.strip()
        return ""

    def parse_text(self, text, model_id):
        """
        Parses the text from a new-format README file and returns a dictionary of metadata.
        Parameters:
            text (str): The README file content.
            model_id (str): The model identifier.

        Returns:
            dict: Parsed metadata with detailed sections.
        """
        results = {"model_id": model_id}
        lines = text.splitlines()
        # Title: first markdown header (# )
        for line in lines:
            if line.startswith("# "):
                results["title"] = line[2:].strip()
                break
        # Extract sections using expected header markers:
        results["description"] = (
            self._extract_section(text, "**Description**") or "No description provided."
        )
        results["model_incorporation_date"] = self._extract_section(
            text, "**Model Incorporation Date:**"
        )
        results["identifiers"] = self._extract_section(text, "## Identifiers")
        results["domain"] = self._extract_section(text, "## Domain")
        results["input"] = self._extract_section(text, "## Input")
        results["output"] = self._extract_section(text, "## Output")
        results["output_columns"] = self._extract_section(
            text, "**Output Columns (up to 10):**"
        )
        results["source_and_deployment"] = self._extract_section(
            text, "## Source and Deployment"
        )
        results["resource_consumption"] = self._extract_section(
            text, "## Resource Consumption"
        )
        results["references"] = self._extract_section(text, "## References")
        results["license"] = self._extract_section(text, "## License")
        results["about_ersilia"] = self._extract_section(text, "## About Ersilia")
        results["github_url"] = "https://github.com/ersilia-os/{}".format(model_id)
        return results


class ReadmeCardResolver(ErsiliaBase):
    """
    Resolves which README parser to use based on the file content.
    """

    def __init__(self, config_json):
        super().__init__(config_json=config_json)
        self.parser_v0 = ReadmeCardVersion0(config_json)
        self.parser_v1 = ReadmeCardVersion1(config_json)
        # Optionally use a local destination directory if provided in config
        self._dest_dir = config_json.get("dest_dir") if config_json else None

    def _raw_readme_url(self, model_id):
        """
        Constructs the raw GitHub URL for the README file of the model.

        Parameters:
            model_id (str): The model identifier.
        Returns:
            str: The raw URL of the README.
        """
        return (
            "https://raw.githubusercontent.com/ersilia-os/{0}/master/README.md".format(
                model_id
            )
        )

    def get_readme_text(self, model_id):
        """
        Gets the README text for a given model, either locally or from GitHub.

        Parameters:
            model_id (str): The model identifier.

        Returns:
            str or None: The README text, or None if not found.
        """
        # Attempt to read the README file locally if _dest_dir is provided
        if self._dest_dir:
            readme_path = os.path.join(self._dest_dir, model_id, "README.md")
            if os.path.exists(readme_path):
                with open(readme_path, "r") as f:
                    return f.read()
        # Otherwise, fetch the README from GitHub
        r = requests.get(self._raw_readme_url(model_id))
        if r.status_code == 200:
            return r.text
        return None

    def resolve_and_parse(self, model_id):
        """
        Determines the README format and parses the text accordingly.

        Parameters:
            model_id (str): The model identifier.

        Returns:
            dict or None: The parsed README metadata.
        """
        text = self.get_readme_text(model_id)
        if not text:
            return None
        if "## Domain" in text or "**Output Columns (up to 10):**" in text:
            return self.parser_v1.parse_text(text, model_id)
        else:
            return self.parser_v0.parse_text(text, model_id)


class ReadmeCard(ErsiliaBase):
    """
    Updated class to handle the README card of a model.
    Supports both the old and new README formats by leveraging the resolver.
    """

    def __init__(self, config_json):
        super().__init__(config_json=config_json)
        self.resolver = ReadmeCardResolver(config_json)

    def get(self, model_id: str = None, slug: str = None) -> dict:
        """
        Retrieves and parses the README file for a model.
        Automatically determines whether the file is in the old or new format.

        Parameters:
            model_id (str): The model identifier.
            slug (str, optional): The model slug (not used here).

        Returns:
            dict: The parsed metadata extracted from the README file.
        """
        if model_id:
            return self.resolver.resolve_and_parse(model_id)
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
            return card
        else:
            return


class LakeCard(ErsiliaBase):
    """
    Class to handle the lake card of a model.

    The lake in ersilia refers to a result storage platform powered by isaura package to
    store repeated result as a cache and allows user to reuse them. It uses HDF5 explorer to
    explore and retrieve information from HDF5 files.

    Parameters
    ----------
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def get(
        self,
        model_id: str = None,
        slug: str = None,
        as_json: bool = False,
    ) -> dict:
        """
        Get the lake card of a model from the HDF5 explorer.

        Parameters
        ----------
        model_id : str, optional
            The ID of the model.
        slug : str, optional
            The slug of the model.
        as_json : bool, optional
            Whether to return the lake card in JSON format.

        Returns
        -------
        dict
            The lake card of the model.
        """
        if model_id is not None:
            if Hdf5Explorer is None:
                self.logger.debug("No lake found")
                return None
            card = Hdf5Explorer(model_id=model_id).info()
            if as_json:
                return json.dumps(card, indent=4)
            else:
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
        card = mc.get(model_id, slug)
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
