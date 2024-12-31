import sys

from ..default import ERSILIA_MODEL_HUB_URL
from .echo import echo


class ModelNotFound(object):
    """
    A class to handle the scenario when a model is not found.

    Attributes
    ----------
    model : object
        The model object that was not found.

    Methods
    -------
    echo()
        Prints an error message and exits the program.
    """

    def __init__(self, model):
        self.model = model

    def echo(self):
        """
        Prints an error message indicating the model was not found and exits the program.
        """
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
    """
    A class to handle the scenario when a model is not found locally.

    Attributes
    ----------
    model_id : str
        The identifier of the model that was not found locally.

    Methods
    -------
    echo()
        Prints an error message and exits the program.
    """

    def __init__(self, model_id):
        self.model_id = model_id

    def echo(self):
        """
        Prints an error message indicating the model was not found locally and exits the program.
        """
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
