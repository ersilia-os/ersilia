import datetime
import os
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ...core.session import Session


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
        mdl.close()
        echo(":no_entry: Model {0} closed".format(mdl.model_id), fg="green")

        old_file_path = "/ersilia/core/current_session.txt"
        new_file_path = os.path.join(os.path.dirname(old_file_path), datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
        os.rename(old_file_path, new_file_path)