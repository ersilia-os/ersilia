from ... import ErsiliaModel
from ..core.session import Session
from ..utils.session import deregister_model_session
 
def close():
  session = Session(config_json=None)
  model_id = session.current_model_id()
  service_class = session.current_service_class()
  if model_id is None:
    raise RuntimeError("No model was served")
  mdl = ErsiliaModel(model_id, service_class=service_class)
  mdl.close()
  deregister_model_session(model_id)
  return ":no_entry: Model {0} closed".format(mdl.model_id)