import json
import os

try:
    import emoji
except:
    emoji = None
import click

from ... import ErsiliaBase
from ...default import (
    API_SCHEMA_FILE,
    APIS_LIST_FILE,
    CARD_FILE,
    MODEL_SIZE_FILE,
    MODEL_SOURCE_FILE,
    PACKMODE_FILE,
    SERVICE_CLASS_FILE,
)
from ...utils.paths import get_metadata_from_base_dir
from .columns_information import ColumnsInformation


class Information(ErsiliaBase):
    """
    Class to handle the information of a models.

    This class provides methods to get various information about a model,
    such as pack mode, service class, model source, API schema, size, metadata, and card.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.repository_folder = os.path.join(
            self._get_bundle_location(model_id=self.model_id)
        )
        self.dest_folder = os.path.join(self._model_path(model_id=model_id))

    def _get_pack_mode(self):
        pack_mode_file = os.path.join(self.dest_folder, PACKMODE_FILE)
        if os.path.exists(pack_mode_file):
            with open(pack_mode_file, "r") as f:
                return f.read().rstrip()
        else:
            return None

    def _get_service_class(self):
        service_class_file = os.path.join(self.repository_folder, SERVICE_CLASS_FILE)
        if os.path.exists(service_class_file):
            with open(service_class_file, "r") as f:
                return f.read().rstrip()
        else:
            return None

    def _get_model_source(self):
        model_source_file = os.path.join((self.dest_folder), MODEL_SOURCE_FILE)
        if os.path.exists(model_source_file):
            with open(model_source_file) as f:
                return f.read().rstrip()
        else:
            return None

    def _get_api_schema(self):
        api_schema_file = os.path.join(self.dest_folder, API_SCHEMA_FILE)
        if os.path.exists(api_schema_file):
            with open(api_schema_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_size(self):
        size_file = os.path.join(self.dest_folder, MODEL_SIZE_FILE)
        if os.path.exists(size_file):
            with open(size_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_metadata(self):
        try:
            data = get_metadata_from_base_dir(self.dest_folder)
        except FileNotFoundError:
            return None
        return data

    def _get_card(self):
        card_file = os.path.join(self.dest_folder, CARD_FILE)
        if os.path.exists(card_file):
            with open(card_file, "r") as f:
                return json.load(f)
        else:
            return None

    def _get_apis_list(self):
        apis_list_file = os.path.join(self.dest_folder, CARD_FILE)
        if os.path.exists(apis_list_file):
            with open(os.path.join(self.repository_folder, APIS_LIST_FILE), "r") as f:
                return [x.rstrip() for x in f.readlines()]
        else:
            return None

    def _get_columns(self):
        columns_data = {}
        api_names = self._get_apis_list()
        for api_name in api_names:
            ci = ColumnsInformation(
                model_id=self.model_id, api_name=api_name, config_json=self.config_json
            )
            data = ci.load()
            columns_data[api_name] = data
        return columns_data

    def get(self) -> dict:
        """
        Get various information about the model.

        Returns
        -------
        dict
            A dictionary containing several information about the model.
        """
        data = {
            "pack_mode": self._get_pack_mode(),
            "service_class": self._get_service_class(),
            "model_source": self._get_model_source(),
            "apis_list": self._get_apis_list(),
            "api_schema": self._get_api_schema(),
            "size": self._get_size(),
            "metadata": self._get_metadata(),
            "card": self._get_card(),
            "columns": self._get_columns(),
        }
        return data


class InformationDisplayer(ErsiliaBase):
    """
    Class to display the information of a model.

    This class provides methods to display various information about a model,
    such as description, identifiers, code, and Docker information.

    Parameters
    ----------
    info_data : dict
        The information data of the model.
    config_json : dict, optional
        Configuration settings in JSON format.
    """

    def __init__(self, info_data, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.info_data = info_data
        self.logger.debug(self.info_data)

    @staticmethod
    def _echo(text, **styles):
        if emoji is not None:
            text = emoji.emojize(text)
        return click.echo(click.style(text, **styles))

    def _description_info(self):
        color = "blue"
        card = self.info_data["card"]
        text = ":rocket: {0}".format(card["Title"]).rstrip(os.linesep)
        self._echo(text, fg=color, bold=True)
        text = "{0}".format(card["Description"]).rstrip(os.linesep)
        self._echo(text, fg=color)
        text = ""
        self._echo(text)

    def _identifiers_info(self):
        color = "green"
        card = self.info_data["card"]
        text = ":person_tipping_hand: Identifiers"
        self._echo(text, fg=color, bold=True)
        text = "Model identifiers: {0}".format(card["Identifier"]).rstrip(os.linesep)
        self._echo(text, fg=color)
        text = "Slug: {0}".format(card["Slug"]).rstrip(os.linesep)
        self._echo(text, fg=color)
        text = ""
        self._echo(text)

    def _code_info(self):
        color = "red"
        card = self.info_data["card"]
        text = ":nerd_face: Code and parameters"
        self._echo(text, fg=color, bold=True)
        text = "GitHub: https://github.com/ersilia-os/{0}".format(card["Identifier"])
        self._echo(text, fg=color)
        if "S3" in card:
            s = card["S3"]
        else:
            s = "-"
        text = "AWS S3: {0}".format(s)
        self._echo(text, fg=color)
        text = ""
        self._echo(text)

    def _docker_info(self):
        try:
            color = "blue"
            card = self.info_data["card"]
            dockerhub_field = card["DockerHub"]
            docker_architecture = card["Docker Architecture"]
            text = ":whale: Docker"
            self._echo(text, fg=color, bold=True)
            text = "Docker Hub: {0}".format(dockerhub_field)
            self._echo(text, fg=color)
            text = "Architectures: {0}".format(",".join(docker_architecture))
            self._echo(text, fg=color)
            text = ""
            self._echo(text)
        except:
            self.logger.warning("No metadata for Docker slots")

    def echo(self):
        """
        Display the information about the model.
        """
        self._description_info()
        self._identifiers_info()
        self._code_info()
        self._docker_info()
        text = "For more information, please visit https://ersilia.io/model-hub"
        self._echo(text, fg="black")
