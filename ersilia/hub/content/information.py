import json
import os
import platform

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
        Display the information about the model using a rich Panel layout.
        """
        from rich.console import Console
        from rich.panel import Panel
        from rich.table import Table
        from rich.text import Text

        _service_class_labels = {
            "pulled_docker": "DockerHub",
            "docker": "Docker (local)",
            "conda": "Conda",
            "venv": "Virtual environment",
            "system": "System Python",
            "hosted": "Hosted",
            "dummy": "Dummy",
        }

        console = Console()
        card = self.info_data.get("card") or {}
        model_source = self.info_data.get("model_source")
        service_class_raw = self.info_data.get("service_class")
        service_class = _service_class_labels.get(service_class_raw, service_class_raw)

        def fmt(value):
            if isinstance(value, list):
                return ", ".join(str(v) for v in value)
            return str(value) if value is not None else "—"

        def fmt_arch(value):
            if not isinstance(value, list):
                return fmt(value)
            current = self._current_arch()
            matched = [a for a in value if current in a.lower() or a.lower() in current]
            return matched[0] if matched else fmt(value)

        table = Table(show_header=False, box=None, padding=(0, 1), expand=True)
        table.add_column("Field", style="bold cyan", no_wrap=True, min_width=24)
        table.add_column("Value", overflow="fold")

        # Origin section
        table.add_row(Text(" Origin", style="bold magenta on grey15"), "")
        if model_source:
            table.add_row("  Fetched from", fmt(model_source))
        if service_class:
            table.add_row("  Service class", fmt(service_class))
        if "DockerHub" in card:
            table.add_row("  Docker Hub", fmt(card["DockerHub"]))
        if "Docker Architecture" in card:
            table.add_row("  Architecture", fmt_arch(card["Docker Architecture"]))
        identifier = card.get("Identifier", "")
        if identifier:
            table.add_row("  GitHub", f"https://github.com/ersilia-os/{identifier}")
        table.add_row("", "")

        sections = [
            ("Overview", ["Identifier", "Slug", "Status", "Task", "Subtask"]),
            ("Description", ["Title", "Description", "Interpretation"]),
            ("Input / Output", ["Input", "Input Dimension", "Input Shape", "Output", "Output Dimension", "Output Shape", "Output Type", "Output Consistency"]),
            ("Deployment", ["Deployment", "Source", "Source Type", "S3"]),
            ("Publication", ["License", "Contributor", "Publication Type", "Publication Year", "Publication", "Source Code"]),
            ("Sizes", ["Model Size", "Environment Size", "Image Size"]),
        ]

        for section_title, fields in sections:
            table.add_row(Text(f" {section_title}", style="bold magenta on grey15"), "")
            for field in fields:
                if field in card:
                    table.add_row(f"  {field}", fmt(card[field]))
            table.add_row("", "")

        title = card.get("Title", "")
        panel_title = f"[bold]{identifier}[/bold]  ·  {title}" if identifier else title
        console.print(Panel(table, title=panel_title, border_style="cyan"))
