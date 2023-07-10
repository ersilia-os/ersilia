from . import BaseAction


class TemplatePreparer(BaseAction):
    def __init__(self, model_id, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )

    def _create_pack(self):
        pass

    def _create_dockerfile(self):
        pass

    def _create_service(self):
        pass
  
    def prepare(self):
        self._create_pack()
        self._create_dockerfile()
        self._create_service()
