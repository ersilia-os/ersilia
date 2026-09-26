import json
import os
import platform

from ... import ErsiliaBase
from ...default import (
    API_SCHEMA_FILE,
    APIS_LIST_FILE,
    CARD_FILE,
    DOCKER_INFO_FILE,
    MODEL_SIZE_FILE,
    MODEL_SOURCE_FILE,
    PACKMODE_FILE,
    SERVICE_CLASS_FILE,
    SERVICE_CLASS_LABELS,
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

    def _get_docker_tag(self):
        docker_info_file = os.path.join(self.dest_folder, DOCKER_INFO_FILE)
        if os.path.exists(docker_info_file):
            with open(docker_info_file, "r") as f:
                data = json.load(f)
            return data.get("tag")
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
            "docker_tag": self._get_docker_tag(),
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
    def _current_arch():
        machine = platform.machine().lower()
        if machine in ("arm64", "aarch64"):
            return "linux/arm64"
        return "linux/amd64"

    def echo(self):
        """
        Display the information about the model as a panel.
        """
        card = self.info_data.get("card") or {}
        service_class = self.info_data.get("service_class")
        origin = [
            ("Fetched from", self.info_data.get("model_source")),
            ("Service", SERVICE_CLASS_LABELS.get(service_class, service_class)),
            ("DockerHub", card.get("DockerHub")),
            ("Version", self.info_data.get("docker_tag")),
        ]
        arch = card.get("Docker Architecture")
        if isinstance(arch, list):
            current = self._current_arch()
            matched = [a for a in arch if current in a.lower() or a.lower() in current]
            arch = matched[0] if matched else arch
        origin.append(("Architecture", arch))
        if card.get("Identifier"):
            origin.append(
                ("GitHub", "https://github.com/ersilia-os/" + card["Identifier"])
            )
        print_card_panel(
            card,
            origin=[(k, v) for k, v in origin if v],
            skip={"DockerHub", "Docker Architecture"},
        )


CARD_SECTIONS = [
    ("Overview", ["Identifier", "Slug", "Status", "Task", "Subtask"]),
    ("Description", ["Title", "Description", "Interpretation"]),
    (
        "Input / Output",
        [
            "Input",
            "Input Dimension",
            "Input Shape",
            "Output",
            "Output Dimension",
            "Output Shape",
            "Output Type",
            "Output Consistency",
        ],
    ),
    (
        "Deployment",
        [
            "Deployment",
            "Source",
            "Source Type",
            "Docker Architecture",
            "DockerHub",
            "S3",
        ],
    ),
    (
        "Publication",
        [
            "License",
            "Contributor",
            "Publication Type",
            "Publication Year",
            "Publication",
            "Source Code",
        ],
    ),
    ("Sizes", ["Model Size", "Environment Size", "Image Size"]),
]
_SIZE_FIELDS = {"Model Size", "Environment Size", "Image Size"}


def print_card_panel(card, origin=(), skip=()):
    """
    Print a model card as a panel, grouped in sections.

    Used by ``ersilia info`` and ``ersilia catalog --card``.

    Parameters
    ----------
    card : dict
        The model card.
    origin : list of (str, str), optional
        Rows for a first "Origin" section (where the model came from locally).
    skip : set of str, optional
        Card fields not to show (e.g. because they are already in ``origin``).
    """
    from rich.text import Text

    from ...utils.echo import fields_table, link, print_panel

    def fmt(field, value):
        if value is None:
            return "—"
        if isinstance(value, list):
            value = ", ".join(str(v) for v in value)
        value = str(value)
        if field in _SIZE_FIELDS and not any(
            u in value.upper() for u in ("MB", "GB", "KB")
        ):
            value = f"{value} MB"
        if value.startswith(("http://", "https://")):
            return link(value)
        return value

    table = fields_table()

    def section(title, rows):
        if not rows:
            return
        if table.row_count:
            table.add_row("", "")
        table.add_row(Text(title, style="not dim"), "")
        for label, value in rows:
            table.add_row(f"  {label}", value)

    section("Origin", [(label, fmt(label, value)) for label, value in origin])
    for title, fields in CARD_SECTIONS:
        section(
            title,
            [
                (field, fmt(field, card[field]))
                for field in fields
                if field in card and field not in skip
            ],
        )
    identifier = card.get("Identifier", "")
    title = f"{identifier} · {card.get('Title', '')}" if identifier else None
    print_panel(table, title=title)
