from ersilia import ErsiliaModel
from ...core.session import Session
from ...utils.session import deregister_model_session
from .. echo import echo
 
def close(model_id):
  session = Session(config_json=None)
  service_class = session.current_service_class()
  if model_id is None:
    raise RuntimeError("No model was served")
  mdl = ErsiliaModel(model_id, service_class=service_class)
  mdl.close()
  deregister_model_session(model_id)
  echo(":no_entry: Model {0} closed".format(model_id), fg="green")
  return