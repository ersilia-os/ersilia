from . import ersilia_cli
from .. import echo
from ...core.session import Session


def current_cmd():
    @ersilia_cli.command(
        short_help="Get identifier of current model",
        help="Get identifier of currently served model",
    )
    def current():
        session = Session(config_json=None)
        model_id = session.current_model_id()
        if model_id is not None:
            echo("Current model identifier: {0}".format(model_id), fg="blue")
        else:
            echo("No model is current served...", fg="yellow")
