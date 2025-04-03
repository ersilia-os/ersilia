import os
import shutil

from ..... import ErsiliaBase


class _Symlinker(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _dest_bundle_symlink(self):
        self.logger.debug("Creating model symlink bundle > dest")
        model_id = self.model_id
        path = self._model_path(model_id)
        model_path = os.path.join(path, "model")
        if os.path.exists(model_path):
            shutil.rmtree(model_path)
        bundle_dir = self._get_bundle_location(model_id)
        src = os.path.join(bundle_dir, "model")
        self.logger.debug("Creating symlink from {0}".format(src))
        self.logger.debug("Creating symlink to {0}".format(model_path))
        os.symlink(src, model_path, target_is_directory=True)

    def _symlinks(self):
        self._dest_bundle_symlink()


class BasePack(_Symlinker):
    """
    Base class for handling FastAPI model packs.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict, optional
        Configuration settings for the pack.
    """

    def __init__(self, model_id, config_json=None):
        _Symlinker.__init__(self, model_id, config_json)
