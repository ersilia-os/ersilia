from .test_services import IOService, CheckService, RunnerService
from .. import ErsiliaBase
from ..utils.logging import logger


class ModelTester(ErsiliaBase):
    def __init__(
            self, 
            model_id, 
            env, 
            type, 
            level, 
            dir
        ):
        ErsiliaBase.__init__(
            self, 
            config_json=None, 
            credentials_json=None
        )
        self.model_id = model_id
        self.env = env
        self.type = type
        self.level = level
        self.dir = dir
        self.ios = IOService(
            self.logger, 
            self._dest_dir,
            self._model_path, 
            self._get_bundle_location, 
            self._get_bentoml_location, 
            self.model_id,
            self.env,
            self.type
        )
        self.checks = CheckService(
            self.logger, 
            self.model_id, 
            self._dest_dir,
            self.dir,
            self.ios,
            self.type
        )
        self.runner = RunnerService(
            self.model_id,
            self.logger, 
            self.ios, 
            self.checks, 
            self._model_path,
            self.env,
            self.type,
            self.level,
            self.dir
        )

    def run(self, output_file=None):
        self.runner.run(output_file)