import os
import json
import tempfile
import requests
from pyairtable import Table
from ... import ErsiliaBase
from ...utils.terminal import run_command
from ...auth.auth import Auth

try:
    from isaura.core.hdf5 import Hdf5Explorer
except:
    Hdf5Explorer = None

from ...default import (
    AIRTABLE_MODEL_HUB_BASE_ID,
    AIRTABLE_MODEL_HUB_TABLE_NAME,
    ISAURA_DIR,
    ISAURA_FILE_TAG,
    ISAURA_FILE_TAG_LOCAL,
    H5_EXTENSION,
    CARD_FILE,
)

AIRTABLE_MAX_ROWS = 100000
AIRTABLE_PAGE_SIZE = 100


class ReadmeCard(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _raw_readme_url(self, model_id):
        url = "https://raw.githubusercontent.com/ersilia-os/{0}/master/README.md".format(
            model_id
        )
        return url

    def _gh_view(self, model_id):
        tmp_folder = tempfile.mkdtemp(prefix="ersilia-")
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
        title = self._title(lines)
        descr = self._description(lines)
        results = {
            "model_id": model_id,
            "title": self._title(lines),
            "description": self._description(lines),
            "mode": self._mode(lines),
            "github_url": self._model_github_url(model_id),
        }
        return results

    def get(self, model_id):
        return self.parse(model_id)


class AirtableCard(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.api_key = self._get_read_only_airtable_api_key()
        self.base_id = AIRTABLE_MODEL_HUB_BASE_ID
        self.table_name = AIRTABLE_MODEL_HUB_TABLE_NAME
        self.max_rows = AIRTABLE_MAX_ROWS
        self.page_size = AIRTABLE_PAGE_SIZE
        self.table = Table(self.api_key, self.base_id, self.table_name)

    @staticmethod
    def _get_read_only_airtable_api_key():
        url = "https://raw.githubusercontent.com/ersilia-os/ersilia/master/config/read_only_keys.json"
        r = requests.get(url)
        data = r.json()
        return data["AIRTABLE_READONLY_API_KEY"]

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

    def find_card_by_mode(self, mode):
        return self._find_card(mode, "Mode")

    def get(self, model_id):
        return self.find_card_by_model_id(model_id)


class LocalCard(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def get(self, model_id):
        model_path = self._model_path(model_id)
        card_path = os.path.join(model_path, CARD_FILE)
        if os.path.exists(card_path):
            with open(card_path, "r") as f:
                card = json.load(f)
            return card
        else:
            return None


class LakeCard(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def get(self, model_id, as_json=False):
        if Hdf5Explorer is None:
            self.logger.debug("No lake found")
            return None
        card = Hdf5Explorer(model_id=model_id).info()
        if as_json:
            return json.dumps(card, indent=4)
        else:
            return card


class ModelCard(object):
    def __init__(self, config_json=None):
        self.lc = LocalCard(config_json=config_json)
        self.ac = AirtableCard(config_json=config_json)
        self.rc = ReadmeCard(config_json=config_json)

    def _get(self, model_id):
        card = self.lc.get(model_id)
        if card is not None:
            return card
        card = self.ac.get(model_id)
        if card is not None:
            return card
        card = self.rc.get(model_id)
        if card is not None:
            return card

    def get(self, model_id, as_json=False):
        card = self._get(model_id)
        if card is None:
            return
        if as_json:
            return json.dumps(card, indent=4)
        else:
            return card
