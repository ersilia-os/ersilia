from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact


@artifacts([JSONArtifact('model')])
class Service(BentoService):
    """ BentoML Service class

        Serves the model in a bentoml based web app

        Attributes:
            input: a json file containing the input for the prediction

    """
    @api(input=JsonInput())
    def predict(self, input):
        """Makes a prediction for a specific input.

        Args:
            input: json file

        Returns:
            Model prediction as a "string" with confidence value as integer
        """
        return "Hello Ersilia!"
