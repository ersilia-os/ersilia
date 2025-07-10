from ... import ErsiliaModel
from ...core.session import Session
from ...hub.content.information import InformationDisplayer


def info(model_id):
    # Provides infomration about the current model.

    # Displays info using the InformationDisplayer class?
    session = Session(config_json=None)
    service_class = session.current_service_class()
    if model_id is None:
        raise RuntimeError("No model was served")
    mdl = ErsiliaModel(model_id, service_class=service_class)
    info = mdl.info()
    InformationDisplayer(info).echo()
