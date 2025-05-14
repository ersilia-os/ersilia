import os

from ...default import AIRTABLE_MODEL_HUB_VIEW_URL
from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


def _read_default_fields(field):
    root = os.path.dirname(os.path.abspath(__file__))
    filename = field.lower().replace(" ", "_")
    file_path = os.path.join(
        root, "..", "..", "hub", "content", "metadata", filename + ".txt"
    )
    with open(file_path, "r") as f:
        valid_field = f.read().split("\n")
    return valid_field


class BaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model information\n"
        self.hints = "Please check Ersilia AirTable to make sure you are providing the right information. This is the AirTable link: {0}".format(
            AIRTABLE_MODEL_HUB_VIEW_URL
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class IdentifierBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model identifier"
        self.hints = "Ersilia model identifiers are 7 alphanumeric characters. They always start with eos, followed by a digit. The eos identifier coincides with the name of the repository. Check our current AirTable to see correct identifiers: {0}".format(
            AIRTABLE_MODEL_HUB_VIEW_URL
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class SlugBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia slug"
        self.hints = "Slug must be a 5-60 chars lowercase single-word unique identifier. Use '-' for linking words if necessary"
        ErsiliaError.__init__(self, self.message, self.hints)


class StatusBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia status"
        self.hints = "Only one of the following status is allowed: {}".format(
            ", ".join(_read_default_fields("Status"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class TitleBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia title"
        self.hints = "Title must be a 1 sentence (10 to 300 chars)"
        ErsiliaError.__init__(self, self.message, self.hints)


class DescriptionBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia description"
        self.hints = "Description must be longer than 200 characters and different from the title"
        ErsiliaError.__init__(self, self.message, self.hints)


class ModeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia mode"
        self.hints = "Only one of the following modes is allowed: {}".format(
            ", ".join(_read_default_fields("Mode"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class SourceBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong source information"
        self.hints = "Only one of the following sources is allowed: {}".format(
            ", ".join(_read_default_fields("Source"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class SourceTypeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong source type information"
        self.hints = "Only one of the following source types is allowed: {}".format(
            ", ".join(_read_default_fields("Source Types"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class TaskBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model task"
        self.hints = "Only one of these tasks is allowed: {}.".format(
            ", ".join(_read_default_fields("Task"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class SubtaskBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model subtask"
        self.hints = "Only one of these subtasks is allowed: {}".format(
            ", ".join(_read_default_fields("Subtask"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class BiomedicalAreaBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong biomedical area"
        self.hints = "Only these biomedical areas are allowed: {}. Biomedical areas must be in list format".format(
            ", ".join(_read_default_fields("Biomedical Area"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class TargetOrganismBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong target organism"
        self.hints = "Only these target organisms are allowed: {}. Target organisms must be in list format".format(
            ", ".join(_read_default_fields("Target Organism"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class InputBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia input"
        self.hints = "Only inputs allowed: {}. Input must be in list format".format(
            ", ".join(_read_default_fields("Input"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class InputShapeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia input shape"
        self.hints = "Only one of the following shapes is allowed: {}".format(
            ", ".join(_read_default_fields("Input Shape"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class InputDimensionBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong input dimension"
        self.hints = "Dimension should be at least 1"
        ErsiliaError.__init__(self, self.message, self.hints)


class OutputBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia output"
        self.hints = "Only one of the following outputs is allowed: {}".format(
            ", ".join(_read_default_fields("Output"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class OutputTypeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia output type"
        self.hints = "Only output types allowed: {}. More than one output type can be added in list format".format(
            ", ".join(_read_default_fields("Output Type"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class OutputShapeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia output shape"
        self.hints = "Only one of the following output shapes is allowed: {}".format(
            ", ".join(_read_default_fields("Output Shape"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class OutputDimensionBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong output dimension"
        self.hints = "Dimension should be at least 1"
        ErsiliaError.__init__(self, self.message, self.hints)


class OutputConsistencyBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong consistency"
        self.hints = (
            "Only one of the following output consistency is allowed: {}".format(
                ", ".join(_read_default_fields("Output Consistency"))
            )
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class TagBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model tag"
        self.hints = "Tags must be in list format and they must be accepted our team. This means that only tags that are already available in Ersilia are allowed. If you want to include a new tag, please open a pull request (PR) on the 'tag.txt' file from the Ersilia repository."
        ErsiliaError.__init__(self, self.message, self.hints)


class LicenseBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong license"
        self.hints = "Listed licenses are: {}. If the model has a license not in this list, please open a PR on the 'license.txt' file in the Ersilia repository".format(
            ", ".join(_read_default_fields("License"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class GithubBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia GitHub URL"
        self.hints = "The model does not seem to be publicly available in Ersilia's GitHub organization profile (ersilia-os). Make sure that a model identifier has been set."
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerhubBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia DockerHub URL"
        self.hints = "The model does not seem to be publicly available in Ersilia's DockerHub organization profile (ersiliaos). Make sure that a model identifier has been set."
        ErsiliaError.__init__(self, self.message, self.hints)


class DockerArchitectureBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Docker architecture"
        self.hints = "Listed Docker architectures are: {}. If you are considering a Docker architecture that is not in this list, please open a PR on the 'docker_architecture.txt' file in the Ersilia repository".format(
            ", ".join(_read_default_fields("Docker Architecture"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class S3BaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia AWS S3 URL"
        self.hints = "The model does not seem to be publicly available in Ersilia's AWS S3 bucket for zipped models. Make sure that a model identifier has been set."
        ErsiliaError.__init__(self, self.message, self.hints)


class BothIdentifiersBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Both identifiers field error"
        self.hints = "Ersilia model identifier and/or slug have not been set yet"
        ErsiliaError.__init__(self, self.message, self.hints)


class PublicationBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Publication field error"
        self.hints = "Publication must be a valid URL"
        ErsiliaError.__init__(self, self.message, self.hints)


class PublicationTypeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong publication type"
        self.hints = "Only one of the following status is allowed: {}".format(
            ", ".join(_read_default_fields("Status"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class PublicationYearBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong publication year"
        self.hints = "Publication year must be valid"
        ErsiliaError.__init__(self, self.message, self.hints)


class SourceCodeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Source Code field error"
        self.hints = "Source Code must be a valid URL"
        ErsiliaError.__init__(self, self.message, self.hints)


class MemoryGbBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Memory Gb field error"
        self.hints = "Memory Gb field must be specified as an integer indicating GB of memory limit"
        ErsiliaError.__init__(self, self.message, self.hints)


class EnvironmentSizeMbBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Environment Size field error"
        self.hints = "Environment Size field must be specified as a valid numeric value indicating MB of model environment and dependencies"
        ErsiliaError.__init__(self, self.message, self.hints)


class ImageSizeMbBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Image Size field error"
        self.hints = "Image Size field must be specified as a valid numeric value indicating MB of model docker image"
        ErsiliaError.__init__(self, self.message, self.hints)


class ModelSizeMbBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Model Size field error"
        self.hints = "Model Size field must be specified as a valid numeric value indicating MB of model directory"
        ErsiliaError.__init__(self, self.message, self.hints)


class ComputationalPerformanceBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Computational Performance field error"
        self.hints = "Computational Performance field must be specified as a valid numeric value indicating the time (s) needed to run the input"
        ErsiliaError.__init__(self, self.message, self.hints)


class ContributorBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Contributor field error"
        self.hints = "Contributor must be a valid github username"
        ErsiliaError.__init__(self, self.message, self.hints)


class IncorporationDateBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Incorporation Date field error"
        self.hints = "Incorporation Date must be a valid ISO date"
        ErsiliaError.__init__(self, self.message, self.hints)


class InterpretationBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Interpretation field error"
        self.hints = "Interpretation must be a string of 10 to 300 chars"
        ErsiliaError.__init__(self, self.hints)


class DeploymentBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong deployment option"
        self.hints = "Only these deployment options are allowed: {}. Deployment options must be in list format".format(
            ", ".join(_read_default_fields("Deployment"))
        )
        ErsiliaError.__init__(self, self.message, self.hints)
