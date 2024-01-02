import os
import json

try:
    import emoji
except:
    emoji = None
import click

from ... import ErsiliaBase
from ...default import (
    PACKMODE_FILE,
    API_SCHEMA_FILE,
    MODEL_SIZE_FILE,
    METADATA_JSON_FILE,
    CARD_FILE,
    SERVICE_CLASS_FILE,
    APIS_LIST_FILE,
)


class Information(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.repository_folder = os.path.join(
            self._get_bento_location(model_id=self.model_id)
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
        metadata_file = os.path.join(self.dest_folder, METADATA_JSON_FILE)
        if os.path.exists(metadata_file):
            with open(metadata_file, "r") as f:
                return json.load(f)
        else:
            return None

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

    def get(self):
        data = {
            "pack_mode": self._get_pack_mode(),
            "service_class": self._get_service_class(),
            "apis_list": self._get_apis_list(),
            "api_schema": self._get_api_schema(),
            "size": self._get_size(),
            "metadata": self._get_metadata(),
            "card": self._get_card(),
        }
        return data


class InformationDisplayer(ErsiliaBase):
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
        self._description_info()
        self._identifiers_info()
        self._code_info()
        self._docker_info()
        text = "For more information, please visit https://ersilia.io/model-hub"
        self._echo(text, fg="black")
