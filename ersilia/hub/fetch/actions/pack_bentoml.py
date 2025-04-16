import os
import sys

from ....default import PACK_METHOD_BENTOML, PACKMETHOD_FILE, PACKMODE_FILE
from ...bundle.repo import DockerfileFile, ServiceFile
from ..pack.bentoml_pack.mode import AVAILABLE_MODES, PackModeDecision
from ..pack.bentoml_pack.runners import get_runner
from . import BaseAction
from .modify import ModelModifier


class ModelPacker(BaseAction):
    """
    Packs a model using BentoML.

    Parameters
    ----------
    model_id : str
        Identifier of the model to be packed.
    mode : str
        Packing mode to be used.
    config_json : dict
        Configuration settings for the packer.

    Methods
    -------
    pack()
        Packs the model using BentoML.
    """

    def __init__(self, model_id: str, mode: str, config_json: dict):
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

    def _register_pack_method(self):
        path = self._get_bundle_location(self.model_id)
        with open(os.path.join(path, PACKMETHOD_FILE), "w") as f:
            self.logger.debug(
                "Writing pack method {0} to file {1}".format(
                    PACK_METHOD_BENTOML, PACKMETHOD_FILE
                )
            )
            f.write(PACK_METHOD_BENTOML)

    def _run(self):
        runner = get_runner(self.pack_mode)(
            model_id=self.model_id, config_json=self.config_json
        )
        runner.run()

    def pack(self):
        """
        Packs the model using BentoML.
        """
        self._setup()
        self._decide_pack_mode()
        self._run()
        self._reset()
        self._modify()
        self._register_pack_method()
