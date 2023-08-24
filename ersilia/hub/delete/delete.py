import os
import shutil
import os.path
from ... import ErsiliaBase
from ...utils.terminal import run_command
from ...utils.environment import Environment
from ...utils.conda import SimpleConda
from ...utils.docker import is_inside_docker
from ..content.catalog import ModelCatalog
from ...db.environments.localdb import EnvironmentDb
from ...db.hubdata.localslugs import SlugDb
from ...db.environments.managers import DockerManager
from ...db.disk.fetched import FetchedModelsManager
from ..bundle.status import ModelStatus

from ...default import ISAURA_FILE_TAG, ISAURA_FILE_TAG_LOCAL


def rmtree(path):
    try:
        shutil.rmtree(path)
    except:
        pass
    try:
        os.unlink(path)
    except:
        pass


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
        self.logger.info("Removing folder {0}".format(folder))
        rmtree(folder)


class ModelLakeDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.path = self._lake_dir

    def delete_if_exists(self, path):
        if os.path.isfile(path):
            os.remove(path)
        if os.path.islink(path):
            os.remove(path)

    def delete_local(self, model_id):
        path = os.path.join(
            self.path, "{0}{1}.h5".format(model_id, ISAURA_FILE_TAG_LOCAL)
        )
        self.logger.debug("Deleting {0}".format(path))
        self.delete_if_exists(path)

    def delete_public(self, model_id):
        path = os.path.join(self.path, "{0}{1}.h5".format(model_id, ISAURA_FILE_TAG))
        self.logger.debug("Deleting {0}".format(path))
        self.delete_if_exists(path)

    def delete(self, model_id):
        self.delete_local(model_id)
        self.delete_public(model_id)


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
        self.logger.info("Removing folder {0}".format(folder))
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
        bento_folder = self._get_bentoml_location(model_id)
        if bento_folder is not None:
            self.logger.info("Removing bento folder first {0}".format(bento_folder))
            rmtree(bento_folder)
            os.makedirs(bento_folder, exist_ok=True)
        self.logger.info("Removing folder {0}".format(folder))
        rmtree(folder)


class ModelBentoDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    @staticmethod
    def _delete_service(service):
        cmd = "echo yes | bentoml delete %s" % service
        run_command(cmd)

    def _delete(self, model_id, keep_latest=True):
        ml = ModelCatalog()
        try:
            catalog = ml.bentoml()
        except:
            self.logger.debug("No BentoML Catalog available")
            catalog = None
        if catalog is None:
            return
        if len(catalog.data) == 0:
            return
        services = [r[1] for r in catalog.data if r[0] == model_id]
        if keep_latest and len(services) > 1:
            services = services[1:]
        for service in services:
            self.logger.info("Removing BentoML service {0}".format(service))
            self._delete_service(service)

    def delete(self, model_id):
        self._delete(model_id, keep_latest=False)

    def clean(self, model_id):
        self._delete(model_id, keep_latest=True)


class ModelSlugDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.slugdb = SlugDb(config_json=config_json)

    def delete(self, model_id):
        self.slugdb.delete_by_model_id(model_id)


class ModelCondaDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.envdb = EnvironmentDb(config_json=config_json)
        self.envdb.table = "conda"

    def _to_delete(self, model_id):
        env = self.envdb.envs_of_model(model_id)
        if (
            len(env) != 1
        ):  # Does not do anything if more than one model depend on the environment.
            return None
        else:
            return list(env)[0]

    def delete(self, model_id):
        envs = self.envdb.envs_of_model(model_id)
        envs = list(set(list(envs) + [model_id]))
        for env in envs:
            models = self.envdb.models_of_env(env)
            models.update([model_id])
            if len(models) == 1:
                self.logger.info("Deleting conda environment {0}".format(env))
                try:
                    conda = SimpleConda()
                    conda.delete(env)
                except:
                    continue
            self.envdb.delete(model_id, env)


class ModelPipDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def pip_uninstall(self, model_id):
        run_command("echo y | pip uninstall %s" % model_id)

    def delete(self, model_id):
        env = Environment()
        if env.has_module(model_id):
            self.logger.info("Uninstalling pip package {0}".format(model_id))
            self.pip_uninstall(model_id)


class ModelDockerDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def delete(self, model_id):
        if is_inside_docker():
            return
        self.logger.info(
            "Removing docker images and stopping containers related to {0}".format(
                model_id
            )
        )
        dm = DockerManager(config_json=self.config_json)
        dm.delete_images(model_id)


class ModelFetchedEntryDeleter(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.fmm = FetchedModelsManager(config_json=config_json)

    def delete(self, model_id):
        self.fmm.delete(model_id)


class TmpCleaner(ErsiliaBase):
    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def delete(self):
        os.rmdir(self._tmp_dir)
        os.makedirs(self._tmp_dir)


class ModelFullDeleter(ErsiliaBase):
    def __init__(self, config_json=None, overwrite=True):
        self.overwrite = overwrite
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)

    def needs_delete(self, model_id):
        ms = ModelStatus().status(model_id)
        for k, v in ms.items():
            if v:
                return True
        if os.path.exists(self._model_path(model_id)):
            return True
        return False

    def delete(self, model_id):
        self.logger.info("Starting delete of model {0}".format(model_id))
        ModelEosDeleter(self.config_json).delete(model_id)
        ModelSlugDeleter(self.config_json).delete(model_id)
        ModelBundleDeleter(self.config_json).delete(model_id)
        ModelBentoDeleter(self.config_json).delete(model_id)
        if self.overwrite:
            ModelCondaDeleter(self.config_json).delete(model_id)
        ModelTmpDeleter(self.config_json).delete(model_id)
        ModelLakeDeleter(self.config_json).delete(model_id)
        ModelPipDeleter(self.config_json).delete(model_id)
        ModelDockerDeleter(self.config_json).delete(model_id)
        ModelFetchedEntryDeleter(self.config_json).delete(model_id)
        self.logger.success("Model {0} deleted successfully".format(model_id))
