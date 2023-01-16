import os
import json
import tempfile
import requests
from ... import ErsiliaBase
from ...utils.terminal import run_command
from ...auth.auth import Auth
from ...db.hubdata.interfaces import AirtableInterface
import validators
from validators import ValidationFailure

from ...utils.exceptions_utils.card_exceptions import (
    SlugBaseInformationError,
    IdentifierBaseInformationError,
    StatusBaseInformationError,
    TitleBaseInformationError,
    DescriptionBaseInformationError,
    ModeBaseInformationError,
    InputBaseInformationError,
    InputShapeBaseInformationError,
    TaskBaseInformationError,
    TagBaseInformationError,
    PublicationBaseInformationError,
    SourceBaseInformationError,
    LicenseBaseInformationError,
    GithubBaseInformationError,
    BothIdentifiersBaseInformationError,
)
from ...utils.identifiers.model import ModelIdentifier

try:
    from isaura.core.hdf5 import Hdf5Explorer
except:
    Hdf5Explorer = None

from ...default import CARD_FILE, METADATA_JSON_FILE


class BaseInformation(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self._identifier = None
        self._slug = None
        self._status = None
        self._title = None
        self._description = None
        self._github = None

    def _is_valid_url(self, url_string: str) -> bool:
        result = validators.url(url_string)
        if isinstance(result, ValidationFailure):
            return False
        return result

    @property
    def identifier(self):
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        mi = ModelIdentifier()
        if not mi.is_valid(new_identifier):
            raise IdentifierBaseInformationError
        self._identifier = new_identifier

    @property
    def slug(self):
        return self._slug

    @slug.setter
    def slug(self, new_slug):
        if new_slug.lower() != new_slug:
            raise SlugBaseInformationError
        if len(new_slug) > 60:
            raise SlugBaseInformationError
        if len(new_slug) < 5:
            raise SlugBaseInformationError
        self._slug = new_slug

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, new_status):
        if new_status not in ["Test", "Ready", "In progress", "To do"]:
            raise StatusBaseInformationError
        self._status = new_status

    @property
    def title(self):
        return self._title

    @title.setter
    def title(self, new_title):
        if len(new_title) > 300:
            raise TitleBaseInformationError
        if len(new_title) < 10:
            raise TitleBaseInformationError
        self._title = new_title

    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, new_description):
        if len(new_description) < 300:
            raise DescriptionBaseInformationError
        if new_description == self._title:
            raise DescriptionBaseInformationError
        self._description = new_description

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, new_mode):
        if new_mode not in ["Pretrained", "Retrained", "In-house", "Online"]:
            raise ModeBaseInformationError
        self._mode = new_mode

    @property
    def input(self):
        return self._input

    @input.setter
    def input(self, new_input):
        if type(new_input) is not list:
            raise InputBaseInformationError
        for inp in new_input:
            if inp not in ["Compound", "Protein", "Text"]:
                raise InputBaseInformationError
        self._input = new_input

    @property
    def input_shape(self):
        return self._input_shape

    @input_shape.setter
    def input_shape(self, new_input_shape):
        if new_input_shape not in [
            "Single",
            "Pair",
            "List",
            "Pair of Lists",
            "List of Lists",
        ]:
            raise InputShapeBaseInformationError
        self._input_shape = new_input_shape

    @property
    def task(self):
        return self._task

    @task.setter
    def task(self, new_task):
        if type(new_task) is not list:
            raise TaskBaseInformationError
        for nt in new_task:
            if nt not in [
                "Classification",
                "Regression",
                "Generative",
                "Embedding",
                "Similarity",
                "Clustering",
                "Dimensionality reduction",
            ]:
                raise TaskBaseInformationError
        self._task = new_task

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, new_output):
        self._output = new_output

    @property
    def interpretation(self):
        return self._interpretation

    @interpretation.setter
    def interpretation(self, new_interpretation):
        self._interpretation = new_interpretation

    @property
    def tag(self):
        return self._tag

    @tag.setter
    def tag(self, new_tag):
        if type(new_tag) is not list:
            raise TagBaseInformationError
        self._tag = new_tag

    @property
    def publication(self):
        return self._publication

    @publication.setter
    def publication(self, new_publication):
        if not self._is_valid_url(new_publication):
            raise PublicationBaseInformationError
        self._publication = new_publication

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, new_source):
        if not self._is_valid_url(new_source):
            raise SourceBaseInformationError
        self._source = new_source

    @property
    def license(self):
        return self._license

    @license.setter
    def license(self, new_license):
        if new_license not in [
            "MIT",
            "GPLv3",
            "LGPL",
            "Apache",
            "BSD-2",
            "BSD-3",
            "Mozilla",
            "CC",
            "Proprietary",
            "None",
        ]:
            raise LicenseBaseInformationError
        self._license = new_license

    @property
    def date(self):
        return self._date

    @date.setter
    def date(self, new_date):
        self._date = new_date

    @property
    def contributor(self):
        return self._contributor

    @contributor.setter
    def contributor(self, new_contributor):
        self._contributor = new_contributor

    @property
    def github(self):
        model_id = self.identifier
        if model_id is None:
            raise GithubBaseInformationError
        self._github = "https://github.com/ersilia-os/{0}".format(model_id)
        return self._github

    @property
    def both_identifiers(self):
        model_id = self.identifier
        slug = self.slug
        if model_id is None or slug is None:
            raise BothIdentifiersBaseInformationError
        self._both_identifiers = (model_id, slug)
        return self._both_identifiers

    def as_dict(self):
        data = {
            "Identifier": self.identifier,
            "Slug": self.slug,
            "Status": self.status,
            "Title": self.title,
            "Description": self.description,
            "Mode": self.mode,
            "Input": self.input,
            "Input Shape": self.input_shape,
            "Task": self.task,
            "Output": self.output,
            "Interpretation": self.interpretation,
            "Tag": self.tag,
            "Publication": self.publication,
            "Source": self.source,
            "License": self.license,
        }
        return data

    def from_dict(self, data):
        self.identifier = data["Identifier"]
        self.slug = data["Slug"]
        self.status = data["Status"]
        self.title = data["Title"]
        self.description = data["Description"]
        self.mode = data["Mode"]
        self.input = data["Input"]
        self.input_shape = data["Input Shape"]
        self.task = data["Task"]
        self.output = data["Output"]
        self.interpretation = data["Interpretation"]
        self.tag = data["Tag"]
        self.publication = data["Publication"]
        self.source = data["Source"]
        self.license = data["License"]


class RepoMetadataFile(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        self.model_id = model_id
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def _github_url(self, org=None, branch=None):
        if org is None:
            org = "ersilia-os"
        if branch is None:
            branch = "main"
        return "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
            org, self.model_id, branch, METADATA_JSON_FILE
        )

    def get_json_file(self, org=None, branch=None):
        url = self._github_url(org=org, branch=branch)
        self.logger.debug("Reading from {0}".format(url))
        r = requests.get(url)
        if r.status_code == 200:
            text = r.content
            return json.loads(text)
        else:
            return None

    def exists(self, org=None, branch=None):
        if self.get_json_file(org=org, branch=branch) is None:
            return False
        else:
            return True

    def read_information(self, org=None, branch=None, json_path=None):
        if json_path is None:
            data = self.get_json_file(org=org, branch=branch)
        else:
            with open(json_path, "r") as f:
                data = json.load(f)
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

    def _find_record(self):
        data = None
        for records in self.table.iterate(
            page_size=self.page_size, max_records=self.max_rows
        ):
            for record in records:
                if self.model_id == record["fields"]["Identifier"]:
                    data = record["fields"]
        return data

    def read_information(self):
        data = self._find_record()
        bi = BaseInformation(config_json=self.config_json)
        bi.from_dict(data)
        return bi

    def write_information(self, data: BaseInformation):
        d = data.as_dict()
        d["GitHub"] = data.github
        self.table.create(d)


class MetadataCard(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def get(self, model_id):
        dest_dir = self._model_path(model_id=model_id)
        self.logger.debug("Trying to get metadata from: {0}".format(dest_dir))
        metadata_json = os.path.join(dest_dir, METADATA_JSON_FILE)
        if os.path.exists(metadata_json):
            with open(metadata_json, "r") as f:
                data = json.load(f)
            return data
        else:
            return None


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
        self.mc = MetadataCard(config_json=config_json)
        self.rc = ReadmeCard(config_json=config_json)

    def _get(self, model_id):
        card = self.lc.get(model_id)
        if card is not None:
            return card
        card = self.ac.get(model_id)
        if card is not None:
            return card
        card = self.mc.get(model_id)
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
