import os
import shutil

from ..... import ErsiliaBase
from ... import MODEL_INSTALL_COMMANDS_FILE


class _Symlinker(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _dest_bundle_symlink(self):
        # TODO: improve function so that it treats other files.
        #       at the moment it only deals with the model folder
        self.logger.debug("Creating model symlink bundle artifacts > dest")
        model_id = self.model_id
        path = self._model_path(model_id)
        model_path = os.path.join(path, "model")
        if os.path.exists(model_path):
            shutil.rmtree(model_path)
        bundle_dir = self._get_bundle_location(model_id)
        src = os.path.join(bundle_dir, model_id, "artifacts")
        os.symlink(src, model_path, target_is_directory=True)
        # env
        env_path = os.path.join(path, "env")
        if os.path.exists(env_path):
            self.logger.debug("Creating env symlink dest <> bundle")
            trg = os.path.join(bundle_dir, "env")
            self.logger.debug(trg)
            shutil.move(env_path, trg)
            os.symlink(trg, env_path, target_is_directory=True)
        # model_install_commands
        model_install_commands_path = os.path.join(path, MODEL_INSTALL_COMMANDS_FILE)
        if not os.path.exists(model_install_commands_path):
            with open(model_install_commands_path, "w"):
                pass
        trg = os.path.join(bundle_dir, MODEL_INSTALL_COMMANDS_FILE)
        self.logger.debug("Creating model_install_commands.sh symlink dest <> bundle")
        shutil.move(model_install_commands_path, trg)
        os.symlink(trg, model_install_commands_path)

    def _symlinks(self):
        self._dest_bundle_symlink()


class BasePack(_Symlinker):
    """
    Base class for handling model packs.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the pack.
    """

    def __init__(self, model_id, config_json):
        _Symlinker.__init__(self, model_id, config_json)
