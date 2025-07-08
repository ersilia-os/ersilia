import json

import click

from ... import ErsiliaModel
from ...core.session import Session
from ...hub.content.information import InformationDisplayer

<<<<<<< HEAD
def info(model_id):
	# Provides infomration about the current model.
	
	# Displays info using the InformationDisplayer class?
  session = Session(config_json=None)
  service_class = session.current_service_class()
  if model_id is None:
      raise RuntimeError("No model was served")
  mdl = ErsiliaModel(model_id, service_class=service_class) 
  print(mdl.info()) # check that this is a dictionary (maybe simplify if comp.)
=======
def info():
    """
    Provides infomration about the current model.
    
    Displays info using the InformationDisplayer class?
    """
    session = Session(config_json=None)
    model_id = session.current_model_id()
    service_class = session.current_service_class()
    if model_id is None:
        raise RuntimeError("No model was served")
    mdl = ErsiliaModel(model_id, service_class=service_class)
    # info = mdl.info()
    # InformationDisplayer(info).echo()
    return mdl.info() # check that this is a dictionary (maybe simplify if comp.)
>>>>>>> 4173947040b958056feb41f472b80bb3f04fdd3e
