from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact


@artifacts([JSONArtifact("model")])
class Service(BentoService):
    """BentoML Service class

    Serves the model as a bentoml

    Attributes:
        input: dummy input

    """

    @api(input=JsonInput())
    def invert(self, input):
        """Inverts a string.

        Args:
            input: json string

        Returns:
            Inverted string
        """
        return "Inverted!"

    @api(input=JsonInput())
    def shuffle(self, input):
        """Shuffles a string.

        Args:
            input: json string

        Returns:
            Shuffled string
        """
        return "Shuffled!"
