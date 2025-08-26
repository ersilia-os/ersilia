from ... import ErsiliaModel
from ...core.session import Session


def info(model_id):
    """
    Provides information about a specified model.

    This command allows users to get detailed information about a current active session,
    including information about Model Identifiers, Code and Parameters, Docker Hub link and Architectures.

    Args
    -------
    model_id (str): ID of the model.

    Returns
    -------
    function: The info command function to be used by the API.
    str: Confirmation message on success or warning message on failure.

    Raises
    -------
    RuntimeError: If no model was served in the current session.
    """
    session = Session(config_json=None)
    service_class = session.current_service_class()
    if model_id is None:
        raise RuntimeError("No model was served")
    mdl = ErsiliaModel(model_id, service_class=service_class)
    info = mdl.info()
    if "card" in info.keys():
        info = info["card"]
    return info
