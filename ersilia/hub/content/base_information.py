import os
import datetime
import validators

try:
    from validators import ValidationFailure
except ImportError:
    from validators import ValidationError as ValidationFailure


from ...utils.exceptions_utils.base_information_exceptions import (
    SlugBaseInformationError,
    IdentifierBaseInformationError,
    StatusBaseInformationError,
    TitleBaseInformationError,
    DescriptionBaseInformationError,
    ModeBaseInformationError,
    SourceBaseInformationError,
    SourceTypeBaseInformationError,
    InputBaseInformationError,
    InputShapeBaseInformationError,
    OutputBaseInformationError,
    OutputTypeBaseInformationError,
    OutputShapeBaseInformationError,
    OutputDimensionBaseInformationError,
    OutputConsistencyBaseInformationError,
    TaskBaseInformationError,
    SubtaskBaseInformationError,
    BiomedicalAreaBaseInformationError,
    TargetOrganismBaseInformationError,
    TagBaseInformationError,
    PublicationBaseInformationError,
    PublicationTypeBaseInformationError,
    PublicationYearBaseInformationError,
    SourceCodeBaseInformationError,
    LicenseBaseInformationError,
    GithubBaseInformationError,
    DockerhubBaseInformationError,
    DockerArchitectureBaseInformationError,
    S3BaseInformationError,
    BothIdentifiersBaseInformationError,
    MemoryGbBaseInformationError,
)
from ...utils.identifiers.model import ModelIdentifier
from ... import ErsiliaBase


class BaseInformation(ErsiliaBase):
    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
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
        self._output_dimension = None
        self._output_consistency = None
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
    def source(self):
        return self._source

    @source.setter
    def source(self, new_source):
        if new_source not in self._read_default_fields("Source"):
            raise SourceBaseInformationError
        self._source = new_source

    @property
    def source_type(self):
        return self._source_type

    @source_type.setter
    def source_type(self, new_source_type):
        if new_source_type not in self._read_default_fields("Source Type"):
            raise SourceTypeBaseInformationError
        self._source_type = new_source_type

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
        if new_input_shape not in self._read_default_fields("Input Shape"):
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
    def subtask(self):
        return self._subtask

    @subtask.setter
    def subtask(self, new_subtask):
        if type(new_subtask) is str:
            new_subtask = [new_subtask]
        if type(new_subtask) is not list:
            raise SubtaskBaseInformationError
        for nt in new_subtask:
            if nt not in self._read_default_fields("Subtask"):
                raise SubtaskBaseInformationError
        self._subtask = new_subtask

    @property
    def biomedical_area(self):
        return self._biomedical_area

    @biomedical_area.setter
    def biomedical_area(self, new_biomedical_area):
        if type(new_biomedical_area) is str:
            new_biomedical_area = [new_biomedical_area]
        if type(new_biomedical_area) is not list:
            raise BiomedicalAreaBaseInformationError
        for nt in new_biomedical_area:
            if nt not in self._read_default_fields("Biomedical Area"):
                raise BiomedicalAreaBaseInformationError
        self._biomedical_area = new_biomedical_area

    @property
    def target_organism(self):
        return self._target_organism

    @target_organism.setter
    def target_organism(self, new_target_organism):
        if type(new_target_organism) is str:
            new_target_organism = [new_target_organism]
        if type(new_target_organism) is not list:
            raise TargetOrganismBaseInformationError
        for nt in new_target_organism:
            if nt not in self._read_default_fields("Target Organism"):
                raise TargetOrganismBaseInformationError
        self._target_organism = new_target_organism

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
    def output_dimension(self):
        return self._output_dimension

    @output_dimension.setter
    def output_dimension(self, new_output_dimension):
        if type(new_output_dimension) is not int:
            raise OutputDimensionBaseInformationError
        if new_output_dimension < 1:
            raise OutputDimensionBaseInformationError
        self._output_dimension = new_output_dimension

    @property
    def output_consistency(self):
        return self._output_consistency

    @output_consistency.setter
    def output_consistency(self, new_output_consistency):
        default_output_consistency = self._read_default_fields("Output Consistency")
        if new_output_consistency not in default_output_consistency:
            raise OutputConsistencyBaseInformationError
        self._output_consistency = new_output_consistency

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
    def publication_type(self):
        return self._publication_type

    @publication_type.setter
    def publication_type(self, new_publication_type):
        if new_publication_type not in self._read_default_fields("Publication Type"):
            raise PublicationTypeBaseInformationError
        self._publication_type = new_publication_type

    @property
    def publication_year(self):
        return self._publication_year

    @publication_year.setter
    def publication_year(self, new_publication_year):
        if type(new_publication_year) is not int:
            raise PublicationYearBaseInformationError
        if new_publication_year < 1900 or new_publication_year > datetime.today("Y"):
            raise PublicationBaseInformationError
        self._publication_year = new_publication_year

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
        self._github = "https://github.com/ersilia-os/{0}".format(model_id)
        return self._github

    @property
    def dockerhub(self):
        return self._dockerhub

    @dockerhub.setter
    def dockerhub(self, new_dockerhub_url):
        if not new_dockerhub_url.startswith("https://hub.docker.com/r/ersiliaos/"):
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
            if d not in self._read_default_fields("Docker Architecture"):
                raise DockerArchitectureBaseInformationError
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

    def as_dict(self):
        data = {
            "Identifier": self.identifier,
            "Slug": self.slug,
            "Status": self.status,
            "Title": self.title,
            "Description": self.description,
            "Mode": self.mode,
            "Source": self.source,
            "Source Type": self.source_type,
            "Input": self.input,
            "Input Shape": self.input_shape,
            "Task": self.task,
            "Subtask": self.subtask,
            "Biomedical Area": self.biomedical_area,
            "Target organism": self.target_organism,
            "Output": self.output,
            "Output Type": self.output_type,
            "Output Shape": self.output_shape,
            "Output Dimension": self.output_dimension,
            "Output Consistency": self.output_consistency,
            "Interpretation": self.interpretation,
            "Tag": self.tag,
            "Publication": self.publication,
            "Publication Type": self.publication_type,
            "Publication Year": self.publication_year,
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

    def _assign(self, var, key, data):
        if key in data:
            var = data[key]
        else:
            var = None

    def from_dict(self, data):
        self._assign(self.identifier, "Identifier", data)
        self._assign(self.slug, "Slug", data)
        self._assign(self.status, "Status", data)
        self._assign(self.title, "Title", data)
        self._assign(self.description, "Description", data)
        self._assign(self.mode, "Mode", data)
        self._assign(self.source, "Source", data)
        self._assign(self.source_type, "Source Type", data)
        self._assign(self.input, "Input", data)
        self._assign(self.input_shape, "Input Shape", data)
        self._assign(self.task, "Task", data)
        self._assign(self.subtask, "Subtask", data)
        self._assign(self.biomedical_area, "Biomedical Area", data)
        self._assign(self.target_organism, "Target Organism", data)
        self._assign(self.output, "Output", data)
        self._assign(self.output_type, "Output Type", data)
        self._assign(self.output_shape, "Output Shape", data)
        self._assign(self.output_dimension, "Output Dimension", data)
        self._assign(self.output_consistency, "Output Consistency", data)
        self._assign(self.interpretation, "Interpretation", data)
        self._assign(self.tag, "Tag", data)
        self._assign(self.publication, "Publication", data)
        self._assign(self.publication_type, "Publication Type", data)
        self._assign(self.publication_year, "Publication Year", data)
        self._assign(self.source_code, "Source Code", data)
        self._assign(self.license, "License", data)
        self._assign(self.contributor, "Contributor", data)
        self._assign(self.dockerhub, "DockerHub", data)
        self._assign(self.docker_architecture, "Docker Architecture", data)
        self._assign(self.s3, "S3", data)
        self._assign(self.memory_gb, "Memory Gb", data)
