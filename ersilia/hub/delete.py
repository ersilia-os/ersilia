from .. import ErsiliaBase
import os
import shutil
import subprocess
from .list import ModelList


class ModelEosDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _model_path(self, model_id):
        folder = os.path.join(self._dest_dir, model_id)
        return folder

    def delete(self, model_id):
        folder = self._model_path(model_id)
        shutil.rmtree(folder)


class ModelBentoDeleter(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    @staticmethod
    def _delete_service(service):
        proc = subprocess.Popen(['echo yes | bentoml delete %s' % service], shell=True)
        proc.wait()

    def delete(self, model_id, keep_latest=True):
        ml = ModelList()
        df = ml.bentoml()
        df = df[df["MODEL_ID"] == model_id]
        if df.shape[0] == 0:
            return
        services = list(df["BENTO_SERVICE"])
        if keep_latest and len(services) > 1:
            services = services[1:]
        for service in services:
            self._delete_service(service)


class TmpCleaner(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)

    def delete(self):
        os.rmdir(self._tmp_dir)
        os.makedirs(self._tmp_dir)
