from ersilia import ErsiliaModel

from ...core.session import Session
from ...utils.session import deregister_model_session
from ..echo import echo


def close(model_id):
    """
    Close the current session of the served model and clean up associated resources.

    Args
    -------
        model_id (str): ID of the model to close.

    Returns
    -------
        bool: True if the model was successfully closed, False otherwise.
    """
    try:
        session = Session(config_json=None)
        service_class = session.current_service_class()
        if model_id is None:
            echo(":person_tipping_hand: No model was served", fg="yellow")
            return False
        mdl = ErsiliaModel(model_id, service_class=service_class)
        mdl.close()
        deregister_model_session(model_id)
        echo(":no_entry: Model {0} closed".format(model_id), fg="green")
        return True
    except Exception as e:
        echo(f":warning: Failed to close model: {e}", fg="red")
        return False
