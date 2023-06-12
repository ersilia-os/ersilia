from ..utils.ports import find_free_port
from ... import ErsiliaBase


# TODO: This is work in progress
class DockerDaemonService(ErsiliaBase):
    def __init__(self, model_id, config_json=None, preferred_port=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        if preferred_port is None:
            self.port = find_free_port()
        else:
            self.port = preferred_port
        self._preferred_port = preferred_port
        self.logger.debug("Staring Docker Daemon service")

    def is_available(self):
        pass

    def serve(self):
        pass

    def api(self):
        pass

    def close(self):
        pass
