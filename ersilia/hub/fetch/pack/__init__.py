import os
import pathlib
import shutil
from .... import ErsiliaBase
from ...bundle.repo import DockerfileFile

from .. import PYTHON_INSTALLS
from ....default import H5_DATA_FILE, ISAURA_FILE_TAG, H5_EXTENSION


class _Symlinker(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _dest_bundle_symlink(self):
        # TODO: improve function so that it treats other files.
        #       at the moment it only deals with the model folder
        # model
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
        # python_installs
        python_installs_path = os.path.join(path, PYTHON_INSTALLS)
        if not os.path.exists(python_installs_path):
            with open(python_installs_path, "w") as f:
                pass
        trg = os.path.join(bundle_dir, PYTHON_INSTALLS)
        self.logger.debug("Creating python_installs.sh symlink dest <> bundle")
        shutil.move(python_installs_path, trg)
        os.symlink(trg, python_installs_path)

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

    def _dest_lake_symlink(self):
        src = os.path.join(self._model_path(self.model_id), H5_DATA_FILE)
        dst = os.path.join(
            self._lake_dir,
            "{0}{1}{2}".format(self.model_id, ISAURA_FILE_TAG, H5_EXTENSION),
        )
        if os.path.exists(src):
            self.logger.debug("Symbolic link from {0}".format(src))
            self.logger.debug("Symbolic link to {0}".format(dst))
            os.symlink(src, dst, target_is_directory=False)

    def _symlinks(self):
        self._bentoml_bundle_symlink()
        self._dest_bundle_symlink()
        self._dest_lake_symlink()


class _Writer(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _write_python_installs(self):
        self.logger.debug("Writing python installs")
        dockerfile = DockerfileFile(self._model_path(self.model_id))
        dis_warn = "--disable-pip-version-check"
        version = dockerfile.get_bentoml_version()
        runs = dockerfile.get_install_commands()["commands"]
        self.logger.debug("Run commands: {0}".format(runs))
        fn = os.path.join(self._model_path(self.model_id), PYTHON_INSTALLS)
        with open(fn, "w") as f:
            for r in runs:
                if r[:3] == "pip":
                    r = r.split(" ")
                    r = " ".join([r[0]] + [dis_warn] + r[1:])
                f.write("{0}{1}".format(r, os.linesep))
            f.write(
                "pip {2} install bentoml=={0}{1}".format(
                    version["version"], os.linesep, dis_warn
                )
            )
        return fn


class BasePack(_Symlinker, _Writer):
    def __init__(self, model_id, config_json):
        _Symlinker.__init__(self, model_id, config_json)
        _Writer.__init__(self, model_id, config_json)
