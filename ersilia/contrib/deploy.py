from ..core.base import ErsiliaBase
from ..app.app import AppBase, StreamlitApp
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

    @staticmethod
    def _app_type(model_id):
        ab = AppBase()
        if ab._is_streamlit(model_id):
            return "streamlit"
        if ab._is_swagger(model_id):
            return "swagger"
        if ab._is_dash(model_id):
            return "dash"

    def _modify_requirements(self, model_id):
        app_type = self._app_type(model_id)
        if app_type == "streamlit":
            R = self._read_bentoml_requirements(model_id)
            r = R[-1].rstrip("\n") + "\n"
            R[-1] = r
            R += ["streamlit==%s\n" % streamlit.__version__]
            requirements = os.path.join(self._get_tmp_directory(model_id), "requirements.txt")
            with open(requirements, "w") as f:
                for r in R:
                    f.write(r)
            return R
        if app_type == "swagger":
            return
        if app_type == "dash":
            return

    def _modify_dockerfile(self, model_id, envport=False):
        app_type = self._app_type(model_id)
        if app_type == "streamlit":
            R = self._read_bentoml_dockerfile(model_id)
            if envport:
                R[-1] = 'CMD streamlit run /bento/app.py --server.port $PORT\n'
            else:
                R[-1] = 'CMD streamlit run /bento/app.py\n'
            dockerfile = os.path.join(self._get_tmp_directory(model_id), "Dockerfile")
            with open(dockerfile, "w") as f:
                for r in R:
                    f.write(r)
            return R
        if app_type == "swagger":
            return # TODO
        if app_type == "dash":
            return # TODO


class Local(DeployBase):

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(self, config_json=config_json, credentials_json=credentials_json)

    def deploy(self, model_id):
        app_type = self._app_type(model_id)
        if app_type == "streamlit":
            app = StreamlitApp()
            app.run(model_id)
        if app_type == "swagger":
            pass
        if app_type == "dash":
            pass


class Heroku(DeployBase):

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(self, config_json=config_json, credentials_json=credentials_json)
        self._login()

    @staticmethod
    def _login():
        subprocess.Popen("heroku container:login", shell=True).wait()

    def _set_tmp(self, model_id):
        self._copy_to_tmp(model_id)
        self._modify_requirements(model_id)
        self._modify_dockerfile(model_id, envport=True)

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
        subprocess.Popen("heroku apps:destroy %s --confirm %s" % (model_id, model_id))


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


class Deployer(object):

    def __init__(self, cloud="heroku", config_json=None, credentials_json=None):
        """Initialize a cloud deployer. For now, only 'heroku' is available."""
        self.cloud = cloud
        self.dep = None
        if cloud == "local":
            self.dep = Local(config_json=config_json, credentials_json=credentials_json)
        if cloud == "heroku":
            self.dep = Heroku(config_json=config_json, credentials_json=credentials_json)
        if cloud == "aws":
            self.dep = Aws(config_json=config_json, credentials_json=credentials_json)
        if cloud == "googlecloud":
            self.dep = GoogleCloud(config_json=config_json, credentials_json=credentials_json)
        if cloud == "azure":
            self.dep = Azure(config_json=config_json, credentials_json=credentials_json)

    def deploy(self, model_id):
        self.dep.deploy(model_id)