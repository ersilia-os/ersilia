from .exceptions import ErsiliaError
from ...default import AIRTABLE_MODEL_HUB_VIEW_URL


class CardErsiliaError(ErsiliaError):
    def __init__(self):
        self.message = "Error occured while running card command"
        self.hints = ""
        super().__init__(self.message, self.hints)


class BaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model information\n"
        self.hints = "Please check Ersilia AirTable to make sure you are providing the right information. This is the AirTable link: {0}".format(
            AIRTABLE_MODEL_HUB_VIEW_URL
        )
        super().__init__(self.message, self.hints)


class IdentifierBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model identifier"
        self.hints = "Ersilia model identifiers are 7 alphanumeric characters. They always start with eos, followed by a digit. Check our current AirTable to see correct identifiers: {0}".format(
            AIRTABLE_MODEL_HUB_VIEW_URL
        )
        super().__init__(self.message, self.hints)


class SlugBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia slug"
        self.hints = "Slug must be a 5-60 chars lowercase single-word unique identifier. Use '-' for linking words if necessary"
        super().__init__(self.message, self.hints)


class StatusBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia status"
        self.hints = "Only status allowed: Test, Ready, In progress, To do"
        super().__init__(self.message, self.hints)


class TitleBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia title"
        self.hints = "Title must be a 1 sentence (10 to 300 chars)"
        super().__init__(self.message, self.hints)


class DescriptionBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia description"
        self.hints = "Description must be longer than 200 characters and different from the title"
        super().__init__(self.message, self.hints)


class ModeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia mode"
        self.hints = "Only modes allowed: Pretrained, Retrained, In-house, Online"
        super().__init__(self.message, self.hints)


class InputBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia input"
        self.hints = (
            "Only inputs allowed: Compound, Protein, Text. Input must be in list format"
        )
        super().__init__(self.message, self.hints)


class InputShapeBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia input shape"
        self.hints = (
            "Only shapes allowed: Single, Pair, List, Pair of Lists, List of Lists"
        )
        super().__init__(self.message, self.hints)


class OutputBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia output"
        self.hints = "Only outputs allowed: Probability, Score, Compound, Descriptor, Vector, Toxicity, IC50"
        super().__init__(self.message, self.hints)


class TaskBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model task"
        self.hints = "Only tasks allowed: Classification, Regression, Generative, Embedding, Similarity, Clustering, Dimensionality reduction. Tasks must be in list format"
        super().__init__(self.message, self.hints)


class TagBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia model tag"
        self.hints = "Tags must be in list format and they must be accepted our team. This means that only tags that are already available in Ersilia are allowed. If you want to include a new tag, please open a pull request (PR) on the 'tag.txt' file from the Ersilia repository."
        super().__init__(self.message, self.hints)


class LicenseBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong license"
        self.hints = "Listed licenses are:MIT, GPLv3, LGPL, Apache, BSD-2, BSD-3, Mozilla, CC, Proprietary, None"
        super().__init__(self.message, self.hints)


class GithubBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Wrong Ersilia GitHub URL"
        self.hints = "The model does not seem to be publicly available in Ersilia's GitHub organization profile (ersilia-os). Make sure that a model identifier has been set."
        super().__init__(self.message, self.hints)


class BothIdentifiersBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Both identifiers field error"
        self.hints = "Ersilia model identifier and/or slug have not been set yet"
        super().__init__(self.message, self.hints)


class PublicationBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Publication field error"
        self.hints = "Publication must be a valid URL"
        super().__init__(self.message, self.hints)


class SourceBaseInformationError(ErsiliaError):
    def __init__(self):
        self.message = "Source field error"
        self.hints = "Source must be a valid URL"
        super().__init__(self.message, self.hints)
