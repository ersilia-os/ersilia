from .echo import echo
from ..default import ERSILIA_MODEL_HUB_URL

import sys


class ModelNotFound(object):
    def __init__(self, model):
        self.model = model

    def echo(self):
        echo(
            "Model not found... {0} is not a valid model identifier".format(
                self.model.text
            ),
            fg="red",
        )
        echo(
            "Find valid identifiers in the Ersilia Model Hub: {0}".format(
                ERSILIA_MODEL_HUB_URL
            )
        )
        sys.exit(0)


class ModelNotInLocal(object):
    def __init__(self, model_id):
        self.model_id = model_id

    def echo(self):
        echo(
            "Model {0} could not be found in local device".format(self.model_id),
            fg="red",
        )
        echo(
            "Please fetch the model from the Ersilia Model Hub: ersilia fetch {0}".format(
                self.model_id
            )
        )
        sys.exit(0)
