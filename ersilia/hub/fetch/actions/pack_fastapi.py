import os

from ....default import PACK_METHOD_FASTAPI, PACKMETHOD_FILE, PACKMODE_FILE
from ..pack.fastapi_pack.mode import AVAILABLE_MODES, PackModeDecision
from ..pack.fastapi_pack.runners import get_runner
from . import BaseAction


class ModelPacker(BaseAction):
    """
    Packs a model using FastAPI (aka ersilia-pack).

    For more information about the ersilia-pack, visit:
    https://github.com/ersilia-os/ersilia-pack.git

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
        Packs the model using FastAPI.
    """

    def __init__(self, model_id: str, mode: str, config_json: dict):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        if mode is not None:
            assert mode in AVAILABLE_MODES
        self.mode = mode

    def _setup(self):
        self.folder = self._model_path(self.model_id)

    def _decide_pack_mode(self):
        if self.mode is None:
            pmd = PackModeDecision(self.model_id, config_json=self.config_json)
            self.pack_mode = pmd.decide()
        else:
            self.pack_mode = self.mode

        with open(os.path.join(self.folder, PACKMODE_FILE), "w") as f:
            f.write(self.pack_mode)

    def _run(self):
        runner = get_runner(self.pack_mode)(
            model_id=self.model_id, config_json=self.config_json
        )
        runner.run()

    def _register_pack_method(self):
        path = self._get_bundle_location(self.model_id)
        with open(os.path.join(path, PACKMETHOD_FILE), "w") as f:
            self.logger.debug(
                "Writing pack method {0} to file {1}".format(
                    PACK_METHOD_FASTAPI, PACKMETHOD_FILE
                )
            )
            f.write(PACK_METHOD_FASTAPI)

    def pack(self):
        """
        Packs the model using FastAPI.
        """
        self._setup()
        self._decide_pack_mode()
        self._run()
        self._register_pack_method()
