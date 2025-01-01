from ... import ErsiliaModel
from ...core.session import Session
from ...utils.session import deregister_model_session
from .. import echo
from . import ersilia_cli


def close_cmd():
    """
    Closes the current session of the served model.

    This command allows users to close the current session and clean up any resources.

    Returns
    -------
    function
        The close command function to be used by the CLI and for testing in the pytest.

    Examples
    --------
    .. code-block:: console

        Close the current session:
        $ ersilia close
    """

    # Example usage: ersilia close {MODEL}
    @ersilia_cli.command(short_help="Close model", help="Close model")
    def close():
        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        if model_id is None:
            echo("No model was served")
            return
        mdl = ErsiliaModel(model_id, service_class=service_class)
        mdl.close()
        deregister_model_session(model_id)
        echo(":no_entry: Model {0} closed".format(mdl.model_id), fg="green")

    return close
