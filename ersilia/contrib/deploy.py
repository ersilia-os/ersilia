from ..core.base import ErsiliaBase
from ..app.app import AppBase
import subprocess
import os
import shutil
import streamlit


class DeployBase(ErsiliaBase):

    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=credentials_json)

    def _get_bentoml_directory(self, model_id):
        path = os.path.join(self._abs_path(self.cfg.LOCAL.BENTOML), "repository", model_id)
        tag = sorted(os.listdir(path))[-1]
        path = os.path.join(path, tag)
        return path

    def _get_tmp_directory(self, model_id):
        path_o = self._get_bentoml_directory(model_id)
        path_d = os.path.join(self._tmp_dir, model_id, os.path.basename(path_o))
        return path_d

    def _read_bentoml_dockerfile(self, model_id):
        dockerfile = os.path.join(self._get_bentoml_directory(model_id), "Dockerfile")
        with open(dockerfile, "r") as f:
            text = f.readlines()
        return text

    def _read_bentoml_requirements(self, model_id):
        requirements = os.path.join(self._get_bentoml_directory(model_id), "requirements.txt")
        with open(requirements, "r") as f:
            text = f.readlines()
        return text

    def _copy_bentoml_to_tmp(self, model_id):
        path_o = self._get_bentoml_directory(model_id)
        path_d = self._get_tmp_directory(model_id)
        os.makedirs(os.path.join(self._tmp_dir, model_id), exist_ok=True)
        if os.path.exists(path_d):
            shutil.rmtree(path_d)
        shutil.copytree(path_o, path_d)
        return path_d

    def _copy_app_to_tmp(self, model_id):
        ab = AppBase()
        app_script = ab.app_script(model_id)
        shutil.copy(app_script, os.path.join(self._get_tmp_directory(model_id), os.path.basename(app_script)))

    def _copy_to_tmp(self, model_id):
        self._copy_bentoml_to_tmp(model_id)
        self._copy_app_to_tmp(model_id)

    def _delete_tmp(self, model_id):
        shutil.rmtree(self._get_tmp_directory(model_id))


class Local(DeployBase):

    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=credentials_json)

    def deploy(self):
        pass


class Heroku(DeployBase):

    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=credentials_json)
        self._login()

    @staticmethod
    def _login():
        subprocess.Popen("heroku container:login", shell=True).wait()

    def _modify_requirements(self, model_id):
        R = self._read_bentoml_requirements(model_id)
        r = R[-1].rstrip("\n") + "\n"
        R[-1] = r
        R += ["streamlit==%s\n" % streamlit.__version__]
        requirements = os.path.join(self._get_tmp_directory(model_id), "requirements.txt")
        with open(requirements, "w") as f:
            for r in R:
                f.write(r)
        return R

    def _modify_dockerfile(self, model_id):
        R = self._read_bentoml_dockerfile(model_id)
        R[-1] = 'CMD streamlit run /bento/app.py --server.port $PORT\n'
        dockerfile = os.path.join(self._get_tmp_directory(model_id), "Dockerfile")
        with open(dockerfile, "w") as f:
            for r in R:
                f.write(r)
        return R

    def _set_tmp(self, model_id):
        self._copy_to_tmp(model_id)
        self._modify_requirements(model_id)
        self._modify_dockerfile(model_id)

    @staticmethod
    def _create_app(model_id):
        subprocess.Popen("heroku create %s" % model_id, shell=True).wait()

    @staticmethod
    def _push(model_id):
        subprocess.Popen("heroku container:push web --app %s" % model_id, shell=True).wait()

    @staticmethod
    def _release(model_id):
        subprocess.Popen("heroku container:release web --app %s" % model_id, shell=True).wait()

    def deploy(self, model_id):
        cwd = os.getcwd()
        self._set_tmp(model_id)
        os.chdir(self._get_tmp_directory(model_id))
        self._create_app(model_id)
        self._push(model_id)
        self._release(model_id)
        os.chdir(cwd)
        self._delete_tmp(model_id)

    @staticmethod
    def destroy(model_id):
        subprocess.Popen("heroku apps:destroy --confirm %s" % model_id)


class Aws(DeployBase):

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(self, config_json=config_json, credentials_json=credentials_json)
        self._login()

    def _login(self):
        pass

    def deploy(self, model_id):
        pass


class GoogleCloud(DeployBase):

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(self, config_json=config_json, credentials_json=credentials_json)
        self._login()

    def _login(self):
        pass

    def deploy(self, model_id):
        pass


class Azure(ErsiliaBase):

    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=credentials_json)
        self._login()

    def _login(self):
        pass

    def deploy(self, model_id):
        pass
