import json
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ...core.session import Session


def info_cmd():
    # Example usage: ersilia info {MODEL}
    @ersilia_cli.command(
        short_help="Get model information", help="Get model information"
    )
    def info():
        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        if model_id is None:
            echo("No model was served")
            return
        mdl = ErsiliaModel(model_id, service_class=service_class)
        info = mdl.info()
        print(json.dumps(info, indent=4))
        return info
