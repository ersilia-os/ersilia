import datetime
import os
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ...core.session import Session
from ...core.tracking import get_persistent_file_path, close_persistent_file


def close_cmd():
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
        # Close our persistent tracking file before closing session so we have access to session info for the model through session.json
        if os.path.isfile(get_persistent_file_path()):
            close_persistent_file(mdl.model_id)
        mdl.close()
        echo(":no_entry: Model {0} closed".format(mdl.model_id), fg="green")