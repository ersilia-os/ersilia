import datetime
import os

import validators

try:
    from validators import ValidationFailure
except ImportError:
    from validators import ValidationError as ValidationFailure


from ... import ErsiliaBase
from ...utils.exceptions_utils.base_information_exceptions import (
    BiomedicalAreaBaseInformationError,
    BothIdentifiersBaseInformationError,
    ComputationalPerformanceHundredBaseInformationError,
    ComputationalPerformanceOneBaseInformationError,
    ComputationalPerformanceTenBaseInformationError,
    DescriptionBaseInformationError,
    DockerArchitectureBaseInformationError,
    DockerhubBaseInformationError,
    EnvironmentSizeMbBaseInformationError,
    GithubBaseInformationError,
    IdentifierBaseInformationError,
    ImageSizeMbBaseInformationError,
    InputBaseInformationError,
    InputShapeBaseInformationError,
    LicenseBaseInformationError,
    ModeBaseInformationError,
    OutputBaseInformationError,
    OutputConsistencyBaseInformationError,
    OutputDimensionBaseInformationError,
    OutputShapeBaseInformationError,
    OutputTypeBaseInformationError,
    PublicationBaseInformationError,
    PublicationTypeBaseInformationError,
    PublicationYearBaseInformationError,
    S3BaseInformationError,
    SlugBaseInformationError,
    SourceBaseInformationError,
    SourceCodeBaseInformationError,
    SourceTypeBaseInformationError,
    StatusBaseInformationError,
    SubtaskBaseInformationError,
    TagBaseInformationError,
    TargetOrganismBaseInformationError,
    TaskBaseInformationError,
    TitleBaseInformationError,
)
from ...utils.identifiers.model import ModelIdentifier


class BaseInformation(ErsiliaBase):
    """
    Base class for handling and validating model information.

    Parameters
    ----------
    config_json : dict
        Configuration data in JSON format.
    """

    def __init__(self, config_json=None):
        """
        Initialize the base information object with a provided configuration.

        Parameters
        ----------
        config_json : dict
            A JSON-compatible dictionary containing configuration data.

        Attributes
        ----------
        _github : None
            Placeholder for GitHub URL of the model.
        _identifier : None
            Placeholder for a unique identifier string.
        _slug : None
            Placeholder for a descriptive slug string.
        _status : None
            Placeholder for the current status of the object.
        _title : None
            Placeholder for the model’s title.
        _description : None
            Placeholder for a description of the model.
        _mode : None
            Placeholder for the training mode, one of 'retrained', 'pretrained', 'in-house', or 'online'.
        _task : None
            Placeholder for the primary task associated with the model, such as 'classification', or 'regression'.
        _subtask : None
            Placeholder for the subtask associated with the model, such as 'activity prediction', or 'featurization'.
        _input : None
            Placeholder for input data specifications, such as 'Compound'.
        _input_shape : None
            Placeholder for the shape of the input data.
        _output : None
            Placeholder for output data specifications, such as 'Probability', or 'Compound'.
        _output_type : None
            Placeholder for the type of output data.
        _output_shape : None
            Placeholder for the shape of the output data.
        _output_dimension : None
            Placeholder for dimensional notes about the output.
        _output_consistency : None
            Placeholder for output consistency metrics, one of 'fixed', or 'variable'.
        _interpretation : None
            Placeholder for interpretation details of the model's output.
        _tag : None
            Placeholder for tags associated with the model.
        _biomedical_area : None
            Placeholder for the biomedical area associated with the model, such as 'ADMET'.
        _target_organism : None
            Placeholder for the target organism associated with the model.
        _publication_type : None
            Placeholder for the type of publication associated with the model, one of Preprint or Peer reviewed.
        _publication_year : None
            Placeholder for the year of publication.
        _publication : None
            Placeholder for publication references.
        _source_code : None
            Placeholder for source code metadata.
        _license : None
            Placeholder for license information.
        _contributor : None
            Placeholder for contributor information.
        _dockerhub : None
            Placeholder for Docker Hub repository details.
        _docker_architecture : None
            Placeholder for Docker image architecture details.
        _s3 : None
            Placeholder for related AWS S3 information.
        _source : None
            Placeholder for the source of the model, one of 'Local', or 'Online'
        _source_type: None
            Placeholder for the type of source of the model, one of 'Internal', 'External', or 'Replicated'.
        _environment_size_mb : None
            Placeholder for environment size in megabytes.
        _image_size_mb : None
            Placeholder for Docker image size in megabytes.
        _computational_performance_one : None
            Placeholder for single-run computational performance.
        _computational_performance_ten : None
            Placeholder for ten-run computational performance.
        _computational_performance_hund : None
            Placeholder for hundred-run computational performance.
        """
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self._github = None
        self._identifier = None
        self._slug = None
        self._status = None
        self._title = None
        self._description = None
        self._mode = None
        self._task = None
        self._subtask = None
        self._input = None
        self._input_shape = None
        self._output = None
        self._output_type = None
        self._output_shape = None
        self._output_dimension = None
        self._output_consistency = None
        self._interpretation = None
        self._tag = None
        self._biomedical_area = None
        self._target_organism = None
        self._publication_type = None
        self._publication_year = None
        self._publication = None
        self._source_code = None
        self._license = None
        self._contributor = None
        self._dockerhub = None
        self._docker_architecture = None
        self._s3 = None
        self._source = None
        self._source_type = None
        self._environment_size_mb = None
        self._image_size_mb = None
        self._computational_performance_one = None
        self._computational_performance_ten = None
        self._computational_performance_hund = None

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
        """
        Get the model identifier.

        Returns
        -------
        str
            The model identifier.
        """
        return self._identifier

    @identifier.setter
    def identifier(self, new_identifier):
        """
        Set the model identifier.

        Parameters
        ----------
        new_identifier : str
            The new model identifier.

        Raises
        ------
        IdentifierBaseInformationError
            If the identifier is not valid.
        """
        mi = ModelIdentifier()
        if not mi.is_valid(new_identifier):
            raise IdentifierBaseInformationError
        self._identifier = new_identifier

    @property
    def slug(self):
        """
        Get the model slug.

        Returns
        -------
        str
            The model slug.
        """
        return self._slug

    @slug.setter
    def slug(self, new_slug):
        """
        Set the model slug.

        Parameters
        ----------
        new_slug : str
            The new model slug.

        Raises
        ------
        SlugBaseInformationError
            If the slug is not valid.
        """
        if new_slug.lower() != new_slug:
            raise SlugBaseInformationError
        if len(new_slug) > 60:
            raise SlugBaseInformationError
        if len(new_slug) < 5:
            raise SlugBaseInformationError
        self._slug = new_slug

    @property
    def status(self):
        """
        Get the model status.

        Returns
        -------
        str
            The model status.
        """
        return self._status

    @status.setter
    def status(self, new_status):
        """
        Set the model status.

        Parameters
        ----------
        new_status : str
            The new model status.

        Raises
        ------
        StatusBaseInformationError
            If the status is not valid.
        """
        if new_status not in self._read_default_fields("Status"):
            raise StatusBaseInformationError
        self._status = new_status

    @property
    def title(self):
        """
        Get the model title.

        Returns
        -------
        str
            The model title.
        """
        return self._title

    @title.setter
    def title(self, new_title):
        """
        Set the model title.

        Parameters
        ----------
        new_title : str
            The new model title.

        Raises
        ------
        TitleBaseInformationError
            If the title is not valid.
        """
        if len(new_title) > 300:
            raise TitleBaseInformationError
        if len(new_title) < 10:
            raise TitleBaseInformationError
        self._title = new_title

    @property
    def description(self):
        """
        Get the model description.

        Returns
        -------
        str
            The model description.
        """
        return self._description

    @description.setter
    def description(self, new_description):
        """
        Set the model description.

        Parameters
        ----------
        new_description : str
            The new model description.

        Raises
        ------
        DescriptionBaseInformationError
            If the description is not valid.
        """
        if len(new_description) < 200:
            raise DescriptionBaseInformationError
        if new_description == self._title:
            raise DescriptionBaseInformationError
        self._description = new_description

    @property
    def mode(self):
        """
        Get the model mode.

        Returns
        -------
        str
            The model mode.
        """
        return self._mode

    @mode.setter
    def mode(self, new_mode):
        """
        Set the model mode.

        Parameters
        ----------
        new_mode : str
            The new model mode.

        Raises
        ------
        ModeBaseInformationError
            If the mode is not valid.
        """
        if new_mode not in self._read_default_fields("Mode"):
            raise ModeBaseInformationError
        self._mode = new_mode

    @property
    def source(self):
        """
        Get the model source.

        Returns
        -------
        str
            The model source.
        """
        return self._source

    @source.setter
    def source(self, new_source):
        """
        Set the model source.

        Parameters
        ----------
        new_source : str
            The new model source.

        Raises
        ------
        SourceBaseInformationError
            If the source is not valid.
        """
        if new_source not in self._read_default_fields("Source"):
            raise SourceBaseInformationError
        self._source = new_source

    @property
    def source_type(self):
        """
        Get the model source type.

        Returns
        -------
        str
            The model source type.
        """
        return self._source_type

    @source_type.setter
    def source_type(self, new_source_type):
        """
        Set the model source type.

        Parameters
        ----------
        new_source_type : str
            The new model source type.

        Raises
        ------
        SourceTypeBaseInformationError
            If the source type is not valid.
        """
        if new_source_type not in self._read_default_fields("Source Type"):
            raise SourceTypeBaseInformationError
        self._source_type = new_source_type

    @property
    def input(self):
        """
        Get the model input.

        Returns
        -------
        list
            The model input.
        """
        return self._input

    @input.setter
    def input(self, new_input):
        """
        Set the model input.

        Parameters
        ----------
        new_input : list or str
            The new model input.

        Raises
        ------
        InputBaseInformationError
            If the input is not valid.
        """
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
        """
        Get the model input shape.

        Returns
        -------
        str
            The model input shape.
        """
        return self._input_shape

    @input_shape.setter
    def input_shape(self, new_input_shape):
        """
        Set the model input shape.

        Parameters
        ----------
        new_input_shape : str
            The new model input shape.

        Raises
        ------
        InputShapeBaseInformationError
            If the input shape is not valid.
        """
        if new_input_shape not in self._read_default_fields("Input Shape"):
            raise InputShapeBaseInformationError
        self._input_shape = new_input_shape

    @property
    def task(self):
        """
        Get the model task.

        Returns
        -------
        list
            The model task.
        """
        return self._task

    @task.setter
    def task(self, new_task):
        """
        Set the model task.

        Parameters
        ----------
        new_task : list or str
            The new model task.

        Raises
        ------
        TaskBaseInformationError
            If the task is not valid.
        """
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
        """
        Get the model subtask.

        Returns
        -------
        list
            The model subtask.
        """
        return self._subtask

    @subtask.setter
    def subtask(self, new_subtask):
        """
        Set the model subtask.

        Parameters
        ----------
        new_subtask : list or str
            The new model subtask.

        Raises
        ------
        SubtaskBaseInformationError
            If the subtask is not valid.
        """
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
        """
        Get the model biomedical area.

        Returns
        -------
        list
            The model biomedical area.
        """
        return self._biomedical_area

    @biomedical_area.setter
    def biomedical_area(self, new_biomedical_area):
        """
        Set the model biomedical area.

        Parameters
        ----------
        new_biomedical_area : list or str
            The new model biomedical area.

        Raises
        ------
        BiomedicalAreaBaseInformationError
            If the biomedical area is not valid.
        """
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
        """
        Get the model target organism.

        Returns
        -------
        list
            The model target organism.
        """
        return self._target_organism

    @target_organism.setter
    def target_organism(self, new_target_organism):
        """
        Set the model target organism.

        Parameters
        ----------
        new_target_organism : list or str
            The new model target organism.

        Raises
        ------
        TargetOrganismBaseInformationError
            If the target organism is not valid.
        """
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
        """
        Get the model output.

        Returns
        -------
        list
            The model output.
        """
        return self._output

    @output.setter
    def output(self, new_output):
        """
        Set the model output.

        Parameters
        ----------
        new_output : list or str
            The new model output.

        Raises
        ------
        OutputBaseInformationError
            If the output is not valid.
        """
        if type(new_output) is str:
            new_output = [new_output]
        default_output = self._read_default_fields("Output")
        for no in new_output:
            if no not in default_output:
                raise OutputBaseInformationError
        self._output = new_output

    @property
    def output_type(self):
        """
        Get the model output type.

        Returns
        -------
        list
            The model output type.
        """
        return self._output_type

    @output_type.setter
    def output_type(self, new_output_type):
        """
        Set the model output type.

        Parameters
        ----------
        new_output_type : list or str
            The new model output type.

        Raises
        ------
        OutputTypeBaseInformationError
            If the output type is not valid.
        """
        if type(new_output_type) is str:
            new_output_type = [new_output_type]
        default_output_type = self._read_default_fields("Output Type")
        for no in new_output_type:
            if no not in default_output_type:
                raise OutputTypeBaseInformationError
        self._output_type = new_output_type

    @property
    def output_shape(self):
        """
        Get the model output shape.

        Returns
        -------
        str
            The model output shape.
        """
        return self._output_shape

    @output_shape.setter
    def output_shape(self, new_output_shape):
        """
        Set the model output shape.

        Parameters
        ----------
        new_output_shape : str
            The new model output shape.

        Raises
        ------
        OutputShapeBaseInformationError
            If the output shape is not valid.
        """
        default_output_shape = self._read_default_fields("Output Shape")
        if new_output_shape not in default_output_shape:
            raise OutputShapeBaseInformationError
        self._output_shape = new_output_shape

    @property
    def output_dimension(self):
        """
        Get the model output dimension.

        Returns
        -------
        int
            The model output dimension.
        """
        return self._output_dimension

    @output_dimension.setter
    def output_dimension(self, new_output_dimension):
        """
        Set the model output dimension.

        Parameters
        ----------
        new_output_dimension : int
            The new model output dimension.

        Raises
        ------
        OutputDimensionBaseInformationError
            If the output dimension is not valid.
        """
        if type(new_output_dimension) is not int:
            raise OutputDimensionBaseInformationError
        if new_output_dimension < 1:
            raise OutputDimensionBaseInformationError
        self._output_dimension = new_output_dimension

    @property
    def output_consistency(self):
        """
        Get the model output consistency.

        Returns
        -------
        str
            The model output consistency.
        """
        return self._output_consistency

    @output_consistency.setter
    def output_consistency(self, new_output_consistency):
        """
        Set the model output consistency.

        Parameters
        ----------
        new_output_consistency : str
            The new model output consistency.

        Raises
        ------
        OutputConsistencyBaseInformationError
            If the output consistency is not valid.
        """
        default_output_consistency = self._read_default_fields("Output Consistency")
        if new_output_consistency not in default_output_consistency:
            raise OutputConsistencyBaseInformationError
        self._output_consistency = new_output_consistency

    @property
    def interpretation(self):
        """
        Get the model interpretation.

        Returns
        -------
        str
            The model interpretation.
        """
        return self._interpretation

    @interpretation.setter
    def interpretation(self, new_interpretation):
        """
        Set the model interpretation.

        Parameters
        ----------
        new_interpretation : str
            The new model interpretation.
        """
        self._interpretation = new_interpretation

    @property
    def tag(self):
        """
        Get the model tags.

        Returns
        -------
        list
            The model tags.
        """
        return self._tag

    @tag.setter
    def tag(self, new_tag):
        """
        Set the model tags.

        Parameters
        ----------
        new_tag : list or str
            The new model tags.

        Raises
        ------
        TagBaseInformationError
            If the tags are not valid.
        """
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
        """
        Get the model publication URL.

        Returns
        -------
        str
            The model publication URL.
        """
        return self._publication

    @publication.setter
    def publication(self, new_publication):
        """
        Set the model publication URL.

        Parameters
        ----------
        new_publication : str
            The new model publication URL.

        Raises
        ------
        PublicationBaseInformationError
            If the publication URL is not valid.
        """
        if not self._is_valid_url(new_publication):
            raise PublicationBaseInformationError
        self._publication = new_publication

    @property
    def publication_type(self):
        """
        Get the model publication type.

        Returns
        -------
        str
            The model publication type.
        """
        return self._publication_type

    @publication_type.setter
    def publication_type(self, new_publication_type):
        """
        Set the model publication type.

        Parameters
        ----------
        new_publication_type : str
            The new model publication type.

        Raises
        ------
        PublicationTypeBaseInformationError
            If the publication type is not valid.
        """
        if new_publication_type not in self._read_default_fields("Publication Type"):
            raise PublicationTypeBaseInformationError
        self._publication_type = new_publication_type

    @property
    def publication_year(self):
        """
        Get the model publication year.

        Returns
        -------
        int
            The model publication year.
        """
        return self._publication_year

    @publication_year.setter
    def publication_year(self, new_publication_year):
        """
        Set the model publication year.

        Parameters
        ----------
        new_publication_year : int
            The new model publication year.

        Raises
        ------
        PublicationYearBaseInformationError
            If the publication year is not valid.
        """
        if type(new_publication_year) is not int:
            raise PublicationYearBaseInformationError
        if new_publication_year < 1900 or new_publication_year > datetime.today("Y"):
            raise PublicationBaseInformationError
        self._publication_year = new_publication_year

    @property
    def source_code(self):
        """
        Get the model source code URL.

        Returns
        -------
        str
            The model source code URL.
        """
        return self._source_code

    @source_code.setter
    def source_code(self, new_source_code):
        """
        Set the model source code URL.

        Parameters
        ----------
        new_source_code : str
            The new model source code URL.

        Raises
        ------
        SourceCodeBaseInformationError
            If the source code URL is not valid.
        """
        if not self._is_valid_url(new_source_code):
            raise SourceCodeBaseInformationError
        self._source_code = new_source_code

    @property
    def license(self):
        """
        Get the model license.

        Returns
        -------
        str
            The model license.
        """
        return self._license

    @license.setter
    def license(self, new_license):
        """
        Set the model license.

        Parameters
        ----------
        new_license : str
            The new model license.

        Raises
        ------
        LicenseBaseInformationError
            If the license is not valid.
        """
        if new_license not in self._read_default_fields("License"):
            raise LicenseBaseInformationError
        self._license = new_license

    @property
    def date(self):
        """
        Get the model date.

        Returns
        -------
        str
            The model date.
        """
        return self._date

    @date.setter
    def date(self, new_date):
        """
        Set the model date.

        Parameters
        ----------
        new_date : str
            The new model date.
        """
        self._date = new_date

    @property
    def contributor(self):
        """
        Get the model contributor.

        Returns
        -------
        str
            The model contributor.
        """
        return self._contributor

    @contributor.setter
    def contributor(self, new_contributor):
        """
        Set the model contributor.

        Parameters
        ----------
        new_contributor : str
            The new model contributor.
        """
        self._contributor = new_contributor

    @property
    def github(self):
        """
        Get the model GitHub URL.

        Returns
        -------
        str
            The model GitHub URL.

        Raises
        ------
        GithubBaseInformationError
            If the identifier is not set.
        """
        model_id = self.identifier
        if model_id is None:
            raise GithubBaseInformationError
        self._github = "https://github.com/ersilia-os/{0}".format(model_id)
        return self._github

    @property
    def dockerhub(self):
        """
        Get the model DockerHub URL.

        Returns
        -------
        str
            The model DockerHub URL.
        """
        return self._dockerhub

    @dockerhub.setter
    def dockerhub(self, new_dockerhub_url):
        """
        Set the model DockerHub URL.

        Parameters
        ----------
        new_dockerhub_url : str
            The new model DockerHub URL.

        Raises
        ------
        DockerhubBaseInformationError
            If the DockerHub URL is not valid.
        """
        if not new_dockerhub_url.startswith("https://hub.docker.com/r/ersiliaos/"):
            raise DockerhubBaseInformationError
        self._dockerhub = new_dockerhub_url

    @property
    def docker_architecture(self):
        """
        Get the model Docker architecture.

        Returns
        -------
        list
            The model Docker architecture.
        """
        return self._docker_architecture

    @docker_architecture.setter
    def docker_architecture(self, new_docker_architecture):
        """
        Set the model Docker architecture.

        Parameters
        ----------
        new_docker_architecture : list or str
            The new model Docker architecture.

        Raises
        ------
        DockerArchitectureBaseInformationError
            If the Docker architecture is not valid.
        """
        if type(new_docker_architecture) is str:
            new_docker_architecture = [new_docker_architecture]
        for d in new_docker_architecture:
            if d not in self._read_default_fields("Docker Architecture"):
                raise DockerArchitectureBaseInformationError
        self._docker_architecture = new_docker_architecture

    @property
    def s3(self):
        """
        Get the model S3 URL.

        Returns
        -------
        str
            The model S3 URL.
        """
        return self._s3

    @s3.setter
    def s3(self, new_s3_url):
        """
        Set the model S3 URL.

        Parameters
        ----------
        new_s3_url : str
            The new model S3 URL.

        Raises
        ------
        S3BaseInformationError
            If the S3 URL is not valid.
        """
        if not new_s3_url.startswith(
            "https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/"
        ):
            raise S3BaseInformationError
        self._s3 = new_s3_url

    @property
    def both_identifiers(self):
        """
        Get both the model identifier and slug.

        Returns
        -------
        tuple
            The model identifier and slug.

        Raises
        ------
        BothIdentifiersBaseInformationError
            If either the identifier or slug is not set.
        """
        model_id = self.identifier
        slug = self.slug
        if model_id is None or slug is None:
            raise BothIdentifiersBaseInformationError
        self._both_identifiers = (model_id, slug)
        return self._both_identifiers

    @property
    def environment_size(self):
        """
        Get the model evironment Size in Mb.

        Returns
        -------
        int
            The model evironment Size in Mb.
        """
        return self._environment_size_mb

    @environment_size.setter
    def environment_size(self, new_environment_size):
        """
        Set the environment size in MB.

        Parameters
        ----------
        new_environment_size : int
            The new environment size in MB.

        Raises
        ------
        EnvironmentSizeMbBaseInformationError
            If the environment size value is not valid.
        """
        if not isinstance(new_environment_size, (int, float)):
            raise EnvironmentSizeMbBaseInformationError
        self._environment_size_mb = new_environment_size

    @property
    def image_size_mb(self):
        """Get the image size in megabytes.

        Returns
        -------
        int
            The size of the image in MB.
        """
        return self._image_size_mb

    @image_size_mb.setter
    def image_size_mb(self, new_image_size_mb):
        """Set the image size in megabytes.

        Parameters
        ----------
        new_image_size_mb : int
            The new image size in MB.

        Raises
        ------
        ImageSizeMbBaseInformationError
            If `new_image_size_mb` is not an integer.
        """
        if not isinstance(new_image_size_mb, (int, float)):
            raise ImageSizeMbBaseInformationError
        self._image_size_mb = new_image_size_mb

    @property
    def computational_performance_one(self):
        """Get the computational performance at level one.

        Returns
        -------
        int or float
            The computational performance metric at level one.
        """
        return self._computational_performance_one

    @computational_performance_one.setter
    def computational_performance_one(self, new_value):
        """Set the computational performance at level one.

        Parameters
        ----------
        new_value : int or float
            The new computational performance value.

        Raises
        ------
        ComputationalPerformanceOneBaseInformationError
            If `new_value` is not an int or float.
        """
        if not isinstance(new_value, (int, float)):
            raise ComputationalPerformanceOneBaseInformationError
        self._computational_performance_one = new_value

    @property
    def computational_performance_ten(self):
        """Get the computational performance at level ten.

        Returns
        -------
        int or float
            The computational performance metric at level ten.
        """
        return self._computational_performance_ten

    @computational_performance_ten.setter
    def computational_performance_ten(self, new_value):
        """Set the computational performance at level ten.

        Parameters
        ----------
        new_value : int or float
            The new computational performance value.

        Raises
        ------
        ComputationalPerformanceTenBaseInformationError
            If `new_value` is not an int or float.
        """
        if not isinstance(new_value, (int, float)):
            raise ComputationalPerformanceTenBaseInformationError
        self._computational_performance_ten = new_value

    @property
    def computational_performance_hund(self):
        """Get the computational performance at level hundred.

        Returns
        -------
        int or float
            The computational performance metric at level hundred.
        """
        return self._computational_performance_hund

    @computational_performance_hund.setter
    def computational_performance_hund(self, new_value):
        """Set the computational performance at level hundred.

        Parameters
        ----------
        new_value : int or float
            The new computational performance value.

        Raises
        ------
        ComputationalPerformanceHundredBaseInformationError
            If `new_value` is not an int or float.
        """
        if not isinstance(new_value, (int, float)):
            raise ComputationalPerformanceHundredBaseInformationError
        self._computational_performance_hund = new_value

    def as_dict(self):
        """
        Convert the model information to a dictionary.

        Returns
        -------
        dict
            The model information as a dictionary.
        """
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
            "Environment Size": self.environment_size,
            "Image Size": self.image_size_mb,
            "Computational Performance 1": self.computational_performance_one,
            "Computational Performance 10": self.computational_performance_ten,
            "Computational Performance 100": self.computational_performance_hund,
        }
        data = dict((k, v) for k, v in data.items() if v is not None)
        return data

    def _assign(self, attr_name, key, data):
        setattr(self, attr_name, data[key] if key in data else None)

    def from_dict(self, data):
        """
        Load the model information from a dictionary.

        Parameters
        ----------
        data : dict
            The model information as a dictionary.
        """
        self._assign("identifier", "Identifier", data)
        self._assign("slug", "Slug", data)
        self._assign("status", "Status", data)
        self._assign("title", "Title", data)
        self._assign("description", "Description", data)
        self._assign("mode", "Mode", data)
        self._assign("source", "Source", data)
        self._assign("source_type", "Source Type", data)
        self._assign("input", "Input", data)
        self._assign("input_shape", "Input Shape", data)
        self._assign("task", "Task", data)
        self._assign("subtask", "Subtask", data)
        self._assign("biomedical_area", "Biomedical Area", data)
        self._assign("target_organism", "Target Organism", data)
        self._assign("output", "Output", data)
        self._assign("output_type", "Output Type", data)
        self._assign("output_shape", "Output Shape", data)
        self._assign("output_dimension", "Output Dimension", data)
        self._assign("output_consistency", "Output Consistency", data)
        self._assign("interpretation", "Interpretation", data)
        self._assign("tag", "Tag", data)
        self._assign("publication", "Publication", data)
        self._assign("publication_type", "Publication Type", data)
        self._assign("publication_year", "Publication Year", data)
        self._assign("source_code", "Source Code", data)
        self._assign("license", "License", data)
        self._assign("contributor", "Contributor", data)
        self._assign("dockerhub", "DockerHub", data)
        self._assign("docker_architecture", "Docker Architecture", data)
        self._assign("s3", "S3", data)
        self._assign("environment_size", "Environment Size", data)
        self._assign("image_size_mb", "Image Size", data)
        self._assign(
            "computational_performance_one", "Computational Performance 1", data
        )
        self._assign(
            "computational_performance_ten", "Computational Performance 10", data
        )
        self._assign(
            "computational_performance_hund", "Computational Performance 100", data
        )
