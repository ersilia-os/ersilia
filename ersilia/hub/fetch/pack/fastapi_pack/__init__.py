import os
import shutil

from ..... import ErsiliaBase
from .....default import H5_DATA_FILE, H5_EXTENSION, ISAURA_FILE_TAG


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

    def _dest_lake_symlink(self):
        src = os.path.join(self._model_path(self.model_id), H5_DATA_FILE)
        dst = os.path.join(
            self._lake_dir,
            "{0}{1}{2}".format(self.model_id, ISAURA_FILE_TAG, H5_EXTENSION),
        )
        if os.path.exists(src) and os.path.exists(os.path.dirname(dst)):
            self.logger.debug("Symbolic link from {0}".format(src))
            self.logger.debug("Symbolic link to {0}".format(dst))
            os.symlink(src, dst, target_is_directory=False)
        else:
            self.logger.info(
                "Could not create symbolic link from {0} to {1}".format(src, dst)
            )

    def _symlinks(self):
        self._dest_bundle_symlink()
        self._dest_lake_symlink()


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
