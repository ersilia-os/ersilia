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
            "Model {0} was not found in the Ersilia Model Hub.".format(self.model.text),
            fg="red",
        )
        echo(
            "Check the identifier or slug. Browse the models at {0}".format(
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
        echo("Model {0} is not available locally.".format(self.model_id), fg="red")
        echo("Fetch it first with 'ersilia fetch {0}'.".format(self.model_id))
        sys.exit(0)


# Shared wordings, so a situation reads the same in every command.


def no_model_served(fg="red"):
    """
    Tell the user that no model is served in this terminal.

    Parameters
    ----------
    fg : str, optional
        "red" when the command cannot continue, "yellow" when nothing needed
        to be done (e.g. ``ersilia close``).
    """
    echo("No model is being served in this terminal.", fg=fg)
    echo("Serve one first with 'ersilia serve MODEL'.")


def wrong_extension(allowed, err=False):
    """
    Tell the user that an output file has an unsupported extension.

    Parameters
    ----------
    allowed : list of str
        The accepted extensions, e.g. [".csv", ".h5"].
    err : bool, optional
        Print to stderr.
    """
    if len(allowed) > 1:
        names = ", ".join(allowed[:-1]) + " or " + allowed[-1]
    else:
        names = allowed[0]
    echo("The output file must end in {0}.".format(names), fg="red", err=err)


def report_error(error):
    """
    Print an error the standard way (message, then hint) and exit with code 1.

    Parameters
    ----------
    error : Exception
        The error, typically an ``ErsiliaError``.
    """
    from ..utils.exceptions_utils.throw_ersilia_exception import (
        user_message_and_hints,
    )

    message, hints = user_message_and_hints(error)
    echo(message, fg="red")
    if hints:
        echo(hints)
    sys.exit(1)
