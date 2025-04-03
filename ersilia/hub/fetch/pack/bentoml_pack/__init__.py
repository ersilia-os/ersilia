import os
import shutil

from ..... import ErsiliaBase
from .....default import BENTOML_PATH
from ....bundle.repo import DockerfileFile
from ....delete.delete import ModelBentoDeleter
from ... import MODEL_INSTALL_COMMANDS_FILE


class _Deleter(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.config_json = config_json
        self._delete_bentoml_if_exists()

    def _delete_bentoml_if_exists(self):
        bentoml_path = os.path.join(BENTOML_PATH, "repository", self.model_id)
        if os.path.exists(bentoml_path):
            self.logger.debug(
                "BentoML path exists! Removing it: {0}".format(bentoml_path)
            )
            deleter = ModelBentoDeleter(config_json=self.config_json)
            deleter.delete(model_id=self.model_id)
        self.logger.debug("Trying to remove path: {0}".format(bentoml_path))
        try:
            self.logger.debug("...successfully")
            shutil.rmtree(bentoml_path)
        except:
            self.logger.debug("...but path did not exist!")


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

    def _bentoml_bundle_symlink(self):
        model_id = self.model_id
        src = self._get_bentoml_location(model_id)
        self.logger.debug("BentoML location is {0}".format(src))
        dst_ = os.path.join(self._bundles_dir, model_id)
        os.makedirs(dst_, exist_ok=True)
        dst = os.path.join(dst_, os.path.basename(src))
        shutil.move(src, dst)
        self.logger.debug("Ersilia Bento location is {0}".format(dst))
        self.logger.debug("Building symlinks between {0} and {1}".format(dst, src))
        os.symlink(dst, src, target_is_directory=True)

    def _symlinks(self):
        self._bentoml_bundle_symlink()
        self._dest_bundle_symlink()


class _Writer(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _write_model_install_commands(self):
        self.logger.debug("Writing install commands")
        dockerfile = DockerfileFile(self._model_path(self.model_id))
        dis_warn = "--disable-pip-version-check"
        version = dockerfile.get_bentoml_version()
        runs = dockerfile.get_install_commands()["commands"]
        self.logger.debug("Run commands: {0}".format(runs))
        fn = os.path.join(self._model_path(self.model_id), MODEL_INSTALL_COMMANDS_FILE)
        self.logger.debug("Writing install commands in {0}".format(fn))
        with open(fn, "w") as f:
            for r in runs:
                if r[:3] == "pip":
                    is_pip3 = None
                    if r[:4] == "pip ":
                        r = "python -m pip " + r[4:]
                        is_pip3 = False
                    if r[:4] == "pip3":
                        is_pip3 = True
                    assert is_pip3 is not None
                    r = r.split(" ")
                    if is_pip3:
                        r = " ".join([r[0]] + [dis_warn] + r[1:])
                    else:
                        r = " ".join(r[:3] + [dis_warn] + r[3:])
                f.write("{0}{1}".format(r, os.linesep))
            if version["version"] == "0.11.0":
                cmd = "python -m pip {1} install git+https://github.com/ersilia-os/bentoml-ersilia.git{0}".format(
                    os.linesep, dis_warn
                )
            else:
                cmd = "python -m pip {2} install bentoml=={0}{1}".format(
                    version["version"], os.linesep, dis_warn
                )
            f.write(cmd)
        return fn


class BasePack(_Deleter, _Symlinker, _Writer):
    """
    Base class for handling BentoML model packs.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    config_json : dict
        Configuration settings for the pack.
    """

    def __init__(self, model_id, config_json):
        _Deleter.__init__(self, model_id, config_json)
        _Symlinker.__init__(self, model_id, config_json)
        _Writer.__init__(self, model_id, config_json)
