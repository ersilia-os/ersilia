import os
import sys

from . import BaseAction
from .modify import ModelModifier
from ..pack.mode import PackModeDecision, AVAILABLE_MODES
from ..pack.runners import get_runner
from ...bundle.repo import ServiceFile, DockerfileFile
from ....default import PACKMODE_FILE


class ModelPacker(BaseAction):
    def __init__(self, model_id, mode, config_json):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        if mode is not None:
            assert mode in AVAILABLE_MODES
        self.mode = mode

    def _setup(self):
        folder = self._model_path(self.model_id)
        ServiceFile(folder).rename_service()
        sys.path.insert(0, folder)
        cwd = os.getcwd()
        os.chdir(folder)
        dockerfile = DockerfileFile(folder)
        version = dockerfile.get_bentoml_version()
        if version is None:
            raise Exception
        self.folder = folder
        self.cwd = cwd
        self.dockerfile = dockerfile
        self.version = version

    def _reset(self):
        os.chdir(self.cwd)
        sys.path.remove(self.folder)

    def _decide_pack_mode(self):
        if self.mode is None:
            pmd = PackModeDecision(self.model_id, config_json=self.config_json)
            self.pack_mode = pmd.decide()
        else:
            self.pack_mode = self.mode

        with open(os.path.join(self.folder, PACKMODE_FILE), "w") as f:
            f.write(self.pack_mode)

    def _modify(self):
        mm = ModelModifier(model_id=self.model_id, config_json=self.config_json)
        mm.modify()

    def _run(self):
        runner = get_runner(self.pack_mode)(
            model_id=self.model_id, config_json=self.config_json
        )
        runner.run()

    def pack(self):
        self._setup()
        self._decide_pack_mode()
        self._run()
        self._reset()
        self._modify()
