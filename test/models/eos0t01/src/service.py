from bentoml import BentoService, api, artifacts
from bentoml.adapters import JsonInput
from bentoml.service.artifacts.common import JSONArtifact


@artifacts([JSONArtifact('model')])
class Service(BentoService):

    @api(input=JsonInput())
    def predict(self, input):
        return "Hello Ersilia!"
