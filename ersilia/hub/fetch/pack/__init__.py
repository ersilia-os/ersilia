import os
import pathlib
from .... import ErsiliaBase
from ...bundle.repo import DockerfileFile

from .. import PYTHON_INSTALLS


class _Symlinker(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _dest_bundle_symlink(self):
        #  TODO: Symbolic links between dest folder and bundle or bentoml folders
        #       This will reduce the storage by 1/2.
        #        Additionally, we may want to delete from dest all of the files that are not necessary or
        #       Are outdated.
        pass

    def _bentoml_bundle_symlink(self):
        model_id = self.model_id
        src = self._get_bentoml_location(model_id)
        self.logger.debug("BentoML location is {0}".format(src))
        dst_ = os.path.join(self._bundles_dir, model_id)
        pathlib.Path(dst_).mkdir(parents=True, exist_ok=True)
        dst = os.path.join(dst_, os.path.basename(src))
        self.logger.debug("Destination location is {0}".format(src))
        self.logger.debug("Building symlinks between {0} and {1}".format(src, dst))
        os.symlink(src, dst, target_is_directory=True)


class _Writer(ErsiliaBase):
    def __init__(self, model_id, config_json):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id

    def _write_python_installs(self):
        dockerfile = DockerfileFile(self._model_path(self.model_id))
        dis_warn = "--disable-pip-version-check"
        version = dockerfile.get_bentoml_version()
        runs = dockerfile.get_install_commands()["commands"]
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
