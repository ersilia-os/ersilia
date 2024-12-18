import os
import json
import requests
import yaml
from ... import ErsiliaBase
from ...utils.terminal import run_command
from ...auth.auth import Auth
from ...db.hubdata.interfaces import JsonModelsInterface
import validators

try:
    from validators import ValidationFailure
except ImportError:
    from validators import ValidationError as ValidationFailure

from ...utils.exceptions_utils.card_exceptions import (
    SlugBaseInformationError,
    IdentifierBaseInformationError,
    StatusBaseInformationError,
    TitleBaseInformationError,
    DescriptionBaseInformationError,
    ModeBaseInformationError,
    InputBaseInformationError,
    InputShapeBaseInformationError,
    OutputBaseInformationError,
    OutputTypeBaseInformationError,
    OutputShapeBaseInformationError,
    TaskBaseInformationError,
    TagBaseInformationError,
    PublicationBaseInformationError,
    SourceCodeBaseInformationError,
    LicenseBaseInformationError,
    GithubBaseInformationError,
    DockerhubBaseInformationError,
    DockerArchitectureInformationError,
    S3BaseInformationError,
    BothIdentifiersBaseInformationError,
    MemoryGbBaseInformationError,
)
from ...utils.identifiers.model import ModelIdentifier
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


class BaseInformation(ErsiliaBase):
    """
    Class to handle the base information of a model card.

    A model card contains metadata about a model, such as its identifier, slug, status,
    title, description, mode, input, output, and other relevant information. And This class provides
    methods to validate and set various fields of a model card.

    Parameters
    ----------
    config_json : dict
        Configuration settings in JSON format.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=None
        )
        self._github = None
        self._identifier = None
        self._slug = None
        self._status = None
        self._title = None
        self._description = None
        self._mode = None
        self._task = None
        self._input = None
        self._input_shape = None
        self._output = None
        self._output_type = None
        self._output_shape = None
        self._interpretation = None
        self._tag = None
        self._publication = None
        self._source_code = None
        self._license = None
        self._contributor = None
        self._dockerhub = None
        self._docker_architecture = None
        self._s3 = None
        self._memory_gb = None

    def _is_valid_url(self, url_string: str) -> bool:
        result = validators.url(url_string)
        if isinstance(result, ValidationFailure):
            return False
        return result

    def _read_default_fields(self, field):
        root = os.path.dirname(os.path.abspath(__file__))
        filename = field.lower().replace(" ", "_")
        file_path = os.path.join(root, "metadata", filename + ".txt")
        with open(file_path, "r") as f:
            valid_field = f.read().split("\n")
        return valid_field

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
        if new_status not in self._read_default_fields("Status"):
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
        if len(new_description) < 200:
            raise DescriptionBaseInformationError
        if new_description == self._title:
            raise DescriptionBaseInformationError
        self._description = new_description

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, new_mode):
        if new_mode not in self._read_default_fields("Mode"):
            raise ModeBaseInformationError
        self._mode = new_mode

    @property
    def input(self):
        return self._input

    @input.setter
    def input(self, new_input):
        if type(new_input) is str:
            new_input = [new_input]
        if type(new_input) is not list:
            raise InputBaseInformationError
        for inp in new_input:
            if inp not in self._read_default_fields("Input"):
                raise InputBaseInformationError
        self._input = new_input

    @property
    def input_shape(self):
        return self._input_shape

    @input_shape.setter
    def input_shape(self, new_input_shape):
        if new_input_shape not in self._read_default_fields(
            "Input Shape"
        ):
            raise InputShapeBaseInformationError
        self._input_shape = new_input_shape

    @property
    def task(self):
        return self._task

    @task.setter
    def task(self, new_task):
        if type(new_task) is str:
            new_task = [new_task]
        if type(new_task) is not list:
            raise TaskBaseInformationError
        for nt in new_task:
            if nt not in self._read_default_fields("Task"):
                raise TaskBaseInformationError
        self._task = new_task

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self, new_output):
        if type(new_output) is str:
            new_output = [new_output]
        default_output = self._read_default_fields("Output")
        for no in new_output:
            if no not in default_output:
                raise OutputBaseInformationError
        self._output = new_output

    @property
    def output_type(self):
        return self._output_type

    @output_type.setter
    def output_type(self, new_output_type):
        if type(new_output_type) is str:
            new_output_type = [new_output_type]
        default_output_type = self._read_default_fields("Output Type")
        for no in new_output_type:
            if no not in default_output_type:
                raise OutputTypeBaseInformationError
        self._output_type = new_output_type

    @property
    def output_shape(self):
        return self._output_shape

    @output_shape.setter
    def output_shape(self, new_output_shape):
        default_output_shape = self._read_default_fields("Output Shape")
        if new_output_shape not in default_output_shape:
            raise OutputShapeBaseInformationError
        self._output_shape = new_output_shape

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
        if type(new_tag) is str:
            new_tag = [new_tag]
        if type(new_tag) is not list:
            raise TagBaseInformationError
        default_tags = self._read_default_fields("Tag")
        for nt in new_tag:
            if nt not in default_tags:
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
    def source_code(self):
        return self._source_code

    @source_code.setter
    def source_code(self, new_source_code):
        if not self._is_valid_url(new_source_code):
            raise SourceCodeBaseInformationError
        self._source_code = new_source_code

    @property
    def license(self):
        return self._license

    @license.setter
    def license(self, new_license):
        if new_license not in self._read_default_fields("License"):
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
        self._github = "https://github.com/ersilia-os/{0}".format(
            model_id
        )
        return self._github

    @property
    def dockerhub(self):
        return self._dockerhub

    @dockerhub.setter
    def dockerhub(self, new_dockerhub_url):
        if not new_dockerhub_url.startswith(
            "https://hub.docker.com/r/ersiliaos/"
        ):
            raise DockerhubBaseInformationError
        self._dockerhub = new_dockerhub_url

    @property
    def docker_architecture(self):
        return self._docker_architecture

    @docker_architecture.setter
    def docker_architecture(self, new_docker_architecture):
        if type(new_docker_architecture) is str:
            new_docker_architecture = [new_docker_architecture]
        for d in new_docker_architecture:
            if d not in self._read_default_fields(
                "Docker Architecture"
            ):
                raise DockerArchitectureInformationError
        self._docker_architecture = new_docker_architecture

    @property
    def s3(self):
        return self._s3

    @s3.setter
    def s3(self, new_s3_url):
        if not new_s3_url.startswith(
            "https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/"
        ):
            raise S3BaseInformationError
        self._s3 = new_s3_url

    @property
    def both_identifiers(self):
        model_id = self.identifier
        slug = self.slug
        if model_id is None or slug is None:
            raise BothIdentifiersBaseInformationError
        self._both_identifiers = (model_id, slug)
        return self._both_identifiers

    @property
    def memory_gb(self):
        return self._memory_gb

    @memory_gb.setter
    def memory_gb(self, new_memory_gb):
        if type(new_memory_gb) != int:
            raise MemoryGbBaseInformationError
        self._memory_gb = new_memory_gb

    def as_dict(self) -> dict:
        """
        Convert the base information to a dictionary.

        Returns
        -------
        dict
            The base information as a dictionary.
        """
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
            "Output Type": self.output_type,
            "Output Shape": self.output_shape,
            "Interpretation": self.interpretation,
            "Tag": self.tag,
            "Publication": self.publication,
            "Source Code": self.source_code,
            "License": self.license,
            "Contributor": self.contributor,
            "DockerHub": self.dockerhub,
            "Docker Architecture": self.docker_architecture,
            "S3": self.s3,
            "Memory Gb": self.memory_gb,
        }
        data = dict((k, v) for k, v in data.items() if v is not None)
        return data

    def from_dict(self, data: dict):
        """
        Set the base information from a dictionary.

        Parameters
        ----------
        data : dict
            The dictionary containing the base information.
        """
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
        self.output_type = data["Output Type"]
        self.output_shape = data["Output Shape"]
        self.interpretation = data["Interpretation"]
        self.tag = data["Tag"]
        self.publication = data["Publication"]
        self.source_code = data["Source Code"]
        self.license = data["License"]
        if "Contributor" in data:
            self.contributor = data["Contributor"]
        if "DockerHub" in data:
            self.dockerhub = data["DockerHub"]
        if "Docker Architecture" in data:
            self.docker_architecture = data["Docker Architecture"]
        if "S3" in data:
            self.s3 = data["S3"]
        if "Memory Gb" in data:
            self.memory_gb = data["Memory Gb"]


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
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=None
        )

    def _github_json_url(self, org=None, branch=None):
        if org is None:
            org = "ersilia-os"
        if branch is None:
            branch = "main"
        return (
            "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
                org, self.model_id, branch, METADATA_JSON_FILE
            )
        )

    def _github_yaml_url(self, org=None, branch=None):
        if org is None:
            org = "ersilia-os"
        if branch is None:
            branch = "main"
        return (
            "https://raw.githubusercontent.com/{0}/{1}/{2}/{3}".format(
                org, self.model_id, branch, METADATA_YAML_FILE
            )
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

    def get_json_or_yaml_file(
        self, org: str = None, branch: str = None
    ) -> dict:
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
            self.logger.debug(
                "Trying to get metadata from: {0}".format(dest_dir)
            )
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
        url = "https://raw.githubusercontent.com/ersilia-os/{0}/master/README.md".format(
            model_id
        )
        return url

    def _gh_view(self, model_id):
        tmp_folder = make_temp_dir(prefix="ersilia-")
        tmp_file = os.path.join(tmp_folder, "view.md")
        cmd = "gh repo view {0}/{1} > {2}".format(
            "ersilia-os", model_id, tmp_file
        )
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
        all_models = self.items_all()
        for model in all_models:
            if model["Identifier"] == model_id:
                return model

    def get_card_by_slug(self, slug):
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
