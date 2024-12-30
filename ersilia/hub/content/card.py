import os
import json
import requests
import yaml
from .base_information import BaseInformation
from ... import ErsiliaBase
from ...utils.terminal import run_command
from ...auth.auth import Auth
from ...db.hubdata.interfaces import AirtableInterface
from ...db.hubdata.json_models_interface import JsonModelsInterface
from ...utils.logging import make_temp_dir

try:
    from isaura.core.hdf5 import Hdf5Explorer
except:
    Hdf5Explorer = None

from ...default import (
    CARD_FILE,
    METADATA_JSON_FILE,
    INFORMATION_FILE,
    METADATA_YAML_FILE,
)
from ...utils.paths import get_metadata_from_base_dir


class RepoMetadataFile(ErsiliaBase):
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
            
    def get_json_or_yaml_file(self, org=None, branch=None):
        return self._get_file_content_from_github(org, branch)

    def exists(self, org=None, branch=None):
        if self.get_json_file(org=org, branch=branch) is None:
            return False
        else:
            return True

    def read_information(self, org=None, branch=None, json_or_yaml_path=None):
        if json_or_yaml_path is None:
            data = self.get_json_or_yaml_file(org=org, branch=branch)
        else:
            with open(json_or_yaml_path, "r") as f:
                if json_or_yaml_path.endswith(".json"):
                    data = json.load(f)
                else:
                    data = yaml.safe_load(f)
        bi = BaseInformation(config_json=self.config_json)
        bi.from_dict(data)
        return bi

    def write_information(self, data: BaseInformation, json_path=None):
        data = data.as_dict()
        if json_path is not None:
            with open(json_path, "w") as f:
                json.dump(data, f, indent=4)
        else:
            return data


class AirtableMetadata(AirtableInterface):
    def __init__(self, model_id, config_json=None):
        self.model_id = model_id
        AirtableInterface.__init__(self, config_json=config_json)
        self._warning_empty_row_message = "The AirTable field Identifier was not found! Please check that there are not empty rows."

    def _find_record(self):
        data = None
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                try:
                    if self.model_id == record["fields"]["Identifier"]:
                        data = record["fields"]
                except:
                    self.logger.warning(self._warning_empty_row_message)

        return data

    def _find_airtable_record_id(self):
        rec_id = None
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                try:
                    if self.model_id == record["fields"]["Identifier"]:
                        rec_id = record["id"]
                except:
                    self.logger.warning(self._warning_empty_row_message)
        return rec_id

    def read_information(self):
        data = self._find_record()
        bi = BaseInformation(config_json=self.config_json)
        bi.from_dict(data)
        return bi

    def write_information(self, data: BaseInformation):
        d = data.as_dict()
        d["GitHub"] = data.github
        rec_id = self._find_airtable_record_id()
        if rec_id is None:
            self.logger.debug(
                "Model {0} does not exist in AirTable. Creating new record".format(
                    self.model_id
                )
            )
            self.table.create(d)
        else:
            self.logger.debug(
                "Model {0} exists in AirTable. Updating record".format(self.model_id)
            )
            self.table.update(record_id=rec_id, fields=d)


class ReadmeMetadata(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        self.model_id = model_id
        ErsiliaBase.__init__(self, config_json=config_json)

    def read_information(self):
        self.logger.debug(
            "Cannot read directly from README file. Using AirTable instead"
        )
        am = AirtableMetadata(model_id=self.model_id)
        bi = am.read_information()
        self.logger.info(bi.as_dict())
        return bi

    def write_information(self, data: BaseInformation, readme_path=None):
        d = data.as_dict()
        text = "# {0}\n\n".format(d["Title"])
        text += "{0}\n\n".format(d["Description"].rstrip("\n"))
        text += "## Identifiers\n\n"
        text += "* EOS model ID: `{0}`\n".format(d["Identifier"])
        text += "* Slug: `{0}`\n\n".format(d["Slug"])
        text += "## Characteristics\n\n"
        text += "* Input: `{0}`\n".format(", ".join(d["Input"]))
        text += "* Input Shape: `{0}`\n".format(d["Input Shape"])
        text += "* Task: `{0}`\n".format(", ".join(d["Task"]))
        text += "* Output: `{0}`\n".format(", ".join(d["Output"]))
        text += "* Output Type: `{0}`\n".format(", ".join(d["Output Type"]))
        text += "* Output Shape: `{0}`\n".format(d["Output Shape"])
        text += "* Interpretation: {0}\n\n".format(d["Interpretation"])
        text += "## References\n\n"
        text += "* [Publication]({0})\n".format(d["Publication"])
        text += "* [Source Code]({0})\n".format(d["Source Code"])
        text += "* Ersilia contributor: [{0}](https://github.com/{0})\n\n".format(
            d["Contributor"]
        )
        text += "## Ersilia model URLs\n"
        text += "* [GitHub]({0})\n".format(data.github)
        if "S3" in d:
            text += "* [AWS S3]({0})\n".format(d["S3"])
        if "DockerHub" in d:
            text += "* [DockerHub]({0}) ({1})\n".format(
                d["DockerHub"], ", ".join(d["Docker Architecture"])
            )
        text += "\n"
        text += "## Citation\n\n"
        text += "If you use this model, please cite the [original authors]({0}) of the model and the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff).\n\n".format(
            d["Publication"]
        )
        text += "## License\n\n"
        text += "This package is licensed under a GPL-3.0 license. The model contained within this package is licensed under a {0} license.\n\n".format(
            d["License"]
        )
        text += "Notice: Ersilia grants access to these models 'as is' provided by the original authors, please refer to the original code repository and/or publication if you use the model in your research.\n\n"
        text += "## About Us\n\n"
        text += "The [Ersilia Open Source Initiative](https://ersilia.io) is a Non Profit Organization ([1192266](https://register-of-charities.charitycommission.gov.uk/charity-search/-/charity-details/5170657/full-print)) with the mission is to equip labs, universities and clinics in LMIC with AI/ML tools for infectious disease research.\n\n"
        text += "[Help us](https://www.ersilia.io/donate) achieve our mission!"
        if readme_path is None:
            return text
        else:
            with open(readme_path, "w") as f:
                f.write(text)


class MetadataCard(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def get(self, model_id=None, slug=None):
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

    def _mode(self, lines):
        text = "\n".join(lines)
        return text.split("# ")[1].split("\n")[1].split("#")[0].strip()

    def _model_github_url(self, model_id):
        return "https://github.com/ersilia-os/{0}".format(model_id)

    def parse(self, model_id):
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
            "mode": self._mode(lines),
            "github_url": self._model_github_url(model_id),
        }
        return results

    def get(self, model_id=None, slug=None):
        if model_id:
            return self.parse(model_id)
        else:
            return None


class AirtableCard(AirtableInterface):
    def __init__(self, config_json):
        AirtableInterface.__init__(self, config_json=config_json)

    def _find_card(self, text, field):
        card = None
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                fields = record["fields"]
                if field not in fields:
                    continue
                if text == record["fields"][field]:
                    card = record["fields"]
        return card

    def find_card_by_model_id(self, model_id):
        return self._find_card(model_id, "Identifier")

    def find_card_by_slug(self, slug):
        return self._find_card(slug, "Slug")

    def get(self, model_id=None, slug=None):
        if model_id is not None:
            return self.find_card_by_model_id(model_id)
        elif slug is not None:
            return self.find_card_by_slug(slug)
        else:
            raise ValueError("Either model_id, slug or mode must be provided")


class LocalCard(ErsiliaBase):
    """
    This class provides information on models that have been fetched and are available locally.
    It retrieves and caches information about the models.
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

    def get(self, model_id=None, slug=None):
        """
        This method returns the card for a model. If the model does not exist, it returns None.
        """
        if model_id:
            card = self._load_data(model_id)
            return card
        else:
            return


class LakeCard(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def get(self, model_id=None, slug=None, as_json=False):
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
    def __init__(self, config_json=None):
        JsonModelsInterface.__init__(self, config_json=config_json)

    def get_card_by_model_id(self, model_id):
        all_models = self.items_all()
        for model in all_models:
            if model["Identifier"] == model_id:
                return model
    
    def get_card_by_slug(self, slug):
        all_models = self.items_all()
        for model in all_models:
            if model["Slug"] == slug:
                return model
            
    def get(self, model_id=None, slug=None):
        if model_id is not None:
            return self.get_card_by_model_id(model_id)
        elif slug is not None:
            return self.get_card_by_slug(slug)
        else:
            raise ValueError("Either model_id or slug must be provided")


class ModelCard(object):
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
        ac = AirtableCard(config_json=self.config_json)
        card = ac.get(model_id, slug)
        if card is not None:
            return card
        rc = ReadmeCard(config_json=self.config_json)
        card = rc.get(model_id, slug)
        if card is not None:
            return card

    def get(self, model_id=None, slug=None, as_json=False):
        card = self._get(model_id, slug)
        if card is None:
            return
        if as_json:
            return json.dumps(card, indent=4)
        else:
            return card
    
