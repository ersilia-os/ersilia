from .. import ErsiliaBase


class LakeSynchronizer(ErsiliaBase):

    def __init__(self, model_id, config_json, credentials_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=credentials_json)
        self.model_id = model_id

    def pull(self):
        self.logger.debug("Pulling lake of model {0}".format(self.model_id))
        pass

    def push(self):
        self.logger.debug("Pushing lake of model {0}".format(self.model_id))
        pass

    def publish(self):
        self.pull()
        self.push()
