import json
import click
from . import ersilia_cli
from .. import echo
from ... import ErsiliaModel
from ...core.session import Session
from ...hub.content.information import InformationDisplayer


def info_cmd():
    # Example usage: ersilia info {MODEL}
    @ersilia_cli.command(
        short_help="Get model information", help="Get model information"
    )
    @click.option(
        "--as_json", is_flag=True, default=False, help="Output as JSON format"
    )
    def info(as_json):
        session = Session(config_json=None)
        model_id = session.current_model_id()
        service_class = session.current_service_class()
        if model_id is None:
            echo("No model was served")
            return
        mdl = ErsiliaModel(model_id, service_class=service_class)
        info = mdl.info()
        if as_json:
            echo(json.dumps(info, indent=4))
            return info
        else:
            InformationDisplayer(info).echo()
