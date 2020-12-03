from .. import ErsiliaBase
import os
import shutil
from ..utils.terminal import run_command
from .catalog import ModelCatalog


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
        df = ml.bentoml()
        if df is None:
            return
        df = df[df["MODEL_ID"] == model_id]
        if df.shape[0] == 0:
            return
        services = list(df["BENTO_SERVICE"])
        if keep_latest and len(services) > 1:
            services = services[1:]
        for service in services:
            self._delete_service(service)

    def delete(self, model_id):
        self._delete(model_id, keep_latest=False)

    def clean(self, model_id):
        self._delete(model_id, keep_latest=True)


class ModelPipDeleter(object):

    def __init__(self):
        pass

    def pip_uninstall(self, model_id):
        run_command("echo y | pip uninstall %s" % model_id, quiet=True)

    def delete(self, model_id):
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

    def delete(self, model_id):
        ModelBentoDeleter(self.config_json).delete(model_id)
        ModelEosDeleter(self.config_json).delete(model_id)
        ModelBundleDeleter(self.config_json).delete(model_id)
        ModelTmpDeleter(self.config_json).delete(model_id)
        ModelPipDeleter().delete(model_id)
