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


class VersionNotSupported:
    """
    A class to handle the scenario when the python version used is less than 3.8

    Attributes
    ----------
    versionMajor : int
        The user's major version (e.g 3 in version 3.8).
    versionMinor : int
        The user's minor version (e.g 8 in version 3.8).

    Methods
    -------
    echo()
        prints an error message to the console in red and exits the program.
    """

    def __init__(self, versionMajor, versionMinor):
        self.versionMajor = versionMajor
        self.versionMinor = versionMinor

    def echo(self):
        """
        Prints an error indicating the python version is not supported and exits the program.
        """
        echo(
            "Error: This application requires Python 3.8 or higher.",
            fg="red",
        )
        echo(
            "Current version: {sys.version_info.major}.{sys.version_info.minor}",
            fg="red",
        )
        sys.exit(1)
