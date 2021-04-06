import os
import shutil
from .. import ErsiliaBase
from ..utils.terminal import run_command
from ..utils.environment import Environment
from ..utils.conda import SimpleConda
from .catalog import ModelCatalog
from ..db.environments.localdb import EnvironmentDb
from .status import ModelStatus


class ModelEosDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._dest_dir, model_id)
        return folder

    def delete(self, model_id):
        folder = self._model_path(model_id)
        if not os.path.exists(folder):
            return
        shutil.rmtree(folder)


class ModelTmpDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._tmp_dir, model_id)
        return folder

    def delete(self, model_id):
        folder = self._model_path(model_id)
        if not os.path.exists(folder):
            return
        shutil.rmtree(folder)


class ModelBundleDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._bundles_dir, model_id)
        return folder

    def delete(self, model_id):
        folder = self._model_path(model_id)
        if not os.path.exists(folder):
            return
        shutil.rmtree(folder)


class ModelBentoDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    @staticmethod
    def _delete_service(service):
        cmd = 'echo yes | bentoml delete %s' % service
        run_command(cmd, quiet=True)

    def _delete(self, model_id, keep_latest=True):
        ml = ModelCatalog()
        catalog = ml.bentoml()
        if catalog is None:
            return
        if len(catalog.data) == 0:
            return
        services = [r[1] for r in catalog.data if r[0] == model_id]
        if keep_latest and len(services) > 1:
            services = services[1:]
        for service in services:
            self._delete_service(service)

    def delete(self, model_id):
        self._delete(model_id, keep_latest=False)

    def clean(self, model_id):
        self._delete(model_id, keep_latest=True)


class ModelCondaDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.envdb = EnvironmentDb(config_json=config_json)
        self.envdb.table = "conda"

    def _to_delete(self, model_id):
        env = self.envdb.envs_of_model(model_id)
        if len(env) != 1: # Does not do anything if more than one model depend on the environment.
            return None
        else:
            return list(env)[0]

    def delete(self, model_id):
        envs = self.envdb.envs_of_model(model_id)
        for env in envs:
            models = self.envdb.models_of_env(env)
            if len(models) == 1:
                conda = SimpleConda()
                conda.delete(env)
            self.envdb.delete(model_id, env)


class ModelPipDeleter(object):

    def __init__(self):
        pass

    def pip_uninstall(self, model_id):
        run_command("echo y | pip uninstall %s" % model_id, quiet=True)

    def delete(self, model_id):
        env = Environment()
        if env.has_module(model_id):
            self.pip_uninstall(model_id)


class TmpCleaner(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def delete(self):
        os.rmdir(self._tmp_dir)
        os.makedirs(self._tmp_dir)


class ModelFullDeleter(object):

    def __init__(self, config_json=None):
        self.config_json = config_json

    def needs_delete(self, model_id):
        ms = ModelStatus().status(model_id)
        for k,v in ms.items():
            if v:
                return True
        return False

    def delete(self, model_id):
        ModelEosDeleter(self.config_json).delete(model_id)
        ModelBundleDeleter(self.config_json).delete(model_id)
        ModelBentoDeleter(self.config_json).delete(model_id)
        ModelCondaDeleter(self.config_json).delete(model_id)
        ModelTmpDeleter(self.config_json).delete(model_id)
        ModelPipDeleter().delete(model_id)
