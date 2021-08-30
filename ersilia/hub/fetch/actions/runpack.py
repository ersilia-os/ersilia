import os
import pathlib
from ....utils.terminal import run_command

QUIET = True


class _Symlinker(object):
    def __init__(self):
        pass

    def _dest_bundle_symlink(self):
        #  TODO: Symbolic links between dest folder and bundle or bentoml folders
        #       This will reduce the storage by 1/2.
        #        Additionally, we may want to delete from dest all of the files that are not necessary or
        #       Are outdated.
        pass

    def _bentoml_bundle_symlink(self):
        model_id = self.model_id
        src = self._get_bentoml_location(model_id)
        dst_ = os.path.join(self._bundles_dir, model_id)
        pathlib.Path(dst_).mkdir(parents=True, exist_ok=True)
        dst = os.path.join(dst_, os.path.basename(src))
        os.symlink(src, dst, target_is_directory=True)


class _Writer(object):
    def __init__(self):
        pass

    def _write_python_install(self):
        dis_warn = "--disable-pip-version-check"
        version = dockerfile.get_bentoml_version()
        runs = dockerfile.get_install_commands()["commands"]
        fn = os.path.join(self._model_path(model_id), PYTHON_INSTALLS)
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


class SystemPack(object):
    def __init__(self, model_id):
        self.model_id = model_id

    def setup(self):
        pass

    def run(self):
        folder = self._model_path(self.model_id)
        script_path = os.path.join(folder, self.cfg.HUB.PACK_SCRIPT)
        run_command("python {0}".format(script_path), quiet=QUIET)
        self._bentoml_bundle_symlink(self.model_id)


class VenvPack(object):
    def __init__(self, model_id):
        self.model_id = model_id

    def setup(self):
        pass

    def run(self):
        pass


class CondaPack(object):
    def __init__(self, model_id):
        pass

    def setup(self):
        pass

    def run(self):
        pass


class DockerPack(object):
    def __init__(self, model_id):
        pass

    def setup(self):
        pass

    def run(self):
        pass
