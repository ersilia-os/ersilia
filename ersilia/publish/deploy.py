import os
import shutil
import subprocess

import streamlit

from ..app.app import AppBase, StreamlitApp
from ..core.base import ErsiliaBase


class DeployBase(ErsiliaBase):
    """
    Base class for deployment operations.

    Parameters
    ----------
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.
    """

    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self.docker_org = self.cfg.EXT.DOCKERHUB_ORG

    def _get_bundle_directory(self, model_id):
        path = os.path.join(
            self._bundles_dir, model_id, self._get_latest_bundle_tag(model_id)
        )
        return path

    def _get_tmp_directory(self, model_id):
        path_o = self._get_bundle_directory(model_id)
        path_d = os.path.join(self._tmp_dir, model_id, os.path.basename(path_o))
        return path_d

    def _read_bundle_dockerfile(self, model_id):
        dockerfile = os.path.join(self._get_bundle_directory(model_id), "Dockerfile")
        with open(dockerfile, "r") as f:
            text = f.readlines()
        return text

    def _read_bundle_requirements(self, model_id):
        requirements = os.path.join(
            self._get_bundle_directory(model_id), "requirements.txt"
        )
        with open(requirements, "r") as f:
            text = f.readlines()
        return text

    def _copy_bundle_to_tmp(self, model_id):
        path_o = self._get_bundle_directory(model_id)
        path_d = self._get_tmp_directory(model_id)
        os.makedirs(os.path.join(self._tmp_dir, model_id), exist_ok=True)
        if os.path.exists(path_d):
            shutil.rmtree(path_d)
        shutil.copytree(path_o, path_d)
        return path_d

    def _copy_app_to_tmp(self, model_id):
        ab = AppBase()
        app_script = ab.app_script(model_id)
        shutil.copy(
            app_script,
            os.path.join(
                self._get_tmp_directory(model_id), os.path.basename(app_script)
            ),
        )

    def _copy_to_tmp(self, model_id):
        self._copy_bundle_to_tmp(model_id)
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

    def _was_ersilia_docker_env(self, model_id):
        head = self._read_bundle_dockerfile(model_id)[0].rstrip()
        manifest = "FROM %s" % self.docker_org
        if manifest in head:
            return True
        else:
            return False

    def _modify_requirements(self, model_id):
        if self._was_ersilia_docker_env(model_id):
            return
        app_type = self._app_type(model_id)
        R = self._read_bundle_requirements(model_id)
        r = R[-1].rstrip("\n") + "\n"
        R[-1] = r
        R += ["git+https://github.com/ersilia-os/ersilia\n"]
        if app_type == "streamlit":
            R += ["streamlit==%s\n" % streamlit.__version__]
            requirements = os.path.join(
                self._get_tmp_directory(model_id), "requirements.txt"
            )
            with open(requirements, "w") as f:
                for r in R:
                    f.write(r)
        if app_type == "swagger":
            return
        if app_type == "dash":
            return

    def _modify_dockerfile(self, model_id, envport=False):
        R = self._read_bundle_dockerfile(model_id)
        app_type = self._app_type(model_id)
        if app_type == "streamlit":
            if envport:
                R[-1] = "CMD streamlit run /bento/app.py --server.port $PORT\n"
            else:
                R[-1] = "CMD streamlit run /bento/app.py\n"
            dockerfile = os.path.join(self._get_tmp_directory(model_id), "Dockerfile")
            with open(dockerfile, "w") as f:
                for r in R:
                    f.write(r)
            return R
        if app_type == "swagger":
            return
        if app_type == "dash":
            return


class Local(DeployBase):
    """
    Class for local deployment.

    Parameters
    ----------
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.
    """

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )

    def deploy(self, model_id: str):
        """
        Deploy the model locally.

        Parameters
        ----------
        model_id : str
            The ID of the model to be deployed.
        """
        app_type = self._app_type(model_id)
        if app_type == "streamlit":
            app = StreamlitApp()
            app.run(model_id)
        if app_type == "swagger":
            pass
        if app_type == "dash":
            pass


class Heroku(DeployBase):
    """
    Class for Heroku deployment, in a cloud platform that allows developers to build, run, and operate applications entirely in the cloud.

    Parameters
    ----------
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.

    Examples
    --------
    .. code-block:: python

        deployer = Heroku(
            config_json="path/to/config.json",
            credentials_json="path/to/credentials.json",
        )
        deployer.deploy("model_id")
    """

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self._login()

    @staticmethod
    def _login():
        subprocess.Popen("heroku container:login", shell=True).wait()

    def _set_tmp(self, model_id: str):
        self._copy_to_tmp(model_id)
        self._modify_requirements(model_id)
        self._modify_dockerfile(model_id, envport=True)

    @staticmethod
    def _create_app(model_id: str):
        subprocess.Popen("heroku create %s" % model_id, shell=True).wait()

    @staticmethod
    def _push(model_id: str):
        subprocess.Popen(
            "heroku container:push web --app %s" % model_id, shell=True
        ).wait()

    @staticmethod
    def _release(model_id: str):
        subprocess.Popen(
            "heroku container:release web --app %s" % model_id, shell=True
        ).wait()

    def deploy(self, model_id: str):
        """
        Deploy the model to Heroku.

        This method handles the entire deployment process to Heroku, including setting up temporary directories,
        creating the Heroku app, pushing the Docker container, and releasing the app.

        Parameters
        ----------
        model_id : str
            The ID of the model to be deployed.
        """
        cwd = os.getcwd()
        self._set_tmp(model_id)
        os.chdir(self._get_tmp_directory(model_id))
        self._create_app(model_id)
        self._push(model_id)
        self._release(model_id)
        os.chdir(cwd)
        self._delete_tmp(model_id)

    @staticmethod
    def destroy(model_id: str):
        """
        Destroy the Heroku app.

        This method permanently deletes the Heroku app associated with the given model ID.

        Parameters
        ----------
        model_id : str
            The ID of the model to be destroyed.
        """
        subprocess.Popen("heroku apps:destroy %s --confirm %s" % (model_id, model_id))


class Aws(DeployBase):
    """
    Class for AWS deployment.

    Parameters
    ----------
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.
    """

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self._login()

    def _login(self):
        pass

    def deploy(self, model_id: str):
        """
        Deploy the model to AWS.

        Parameters
        ----------
        model_id : str
            The ID of the model to be deployed.
        """
        pass


class GoogleCloud(DeployBase):
    """
    Class for Google Cloud deployment.

    Parameters
    ----------
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.
    """

    def __init__(self, config_json=None, credentials_json=None):
        DeployBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self._login()

    def _login(self):
        pass

    def deploy(self, model_id: str):
        """
        Deploy the model to Google Cloud.

        Parameters
        ----------
        model_id : str
            The ID of the model to be deployed.
        """
        pass


class Azure(ErsiliaBase):
    """
    Class for Azure deployment.

    Parameters
    ----------
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.
    """

    def __init__(self, config_json=None, credentials_json=None):
        ErsiliaBase.__init__(
            self, config_json=config_json, credentials_json=credentials_json
        )
        self._login()

    def _login(self):
        pass

    def deploy(self, model_id: str):
        """
        Deploy the model to Azure.

        Parameters
        ----------
        model_id : str
            The ID of the model to be deployed.
        """
        pass


class Deployer(object):
    """
    Class for deploying models to various cloud platforms.

    Parameters
    ----------
    cloud : str, optional
        The cloud platform to deploy to. Default is 'heroku'.
    config_json : str, optional
        Path to the configuration JSON file.
    credentials_json : str, optional
        Path to the credentials JSON file.

    Examples
    --------
    .. code-block:: python

        deployer = Deployer(
            cloud="heroku",
            config_json="path/to/config.json",
            credentials_json="path/to/credentials.json",
        )
        deployer.deploy("model_id")
    """

    def __init__(self, cloud="heroku", config_json=None, credentials_json=None):
        """Initialize a cloud deployer. For now, only 'heroku' is available."""
        self.cloud = cloud
        self.dep = None
        if cloud == "local":
            self.dep = Local(config_json=config_json, credentials_json=credentials_json)
        if cloud == "heroku":
            self.dep = Heroku(
                config_json=config_json, credentials_json=credentials_json
            )
        if cloud == "aws":
            self.dep = Aws(config_json=config_json, credentials_json=credentials_json)
        if cloud == "googlecloud":
            self.dep = GoogleCloud(
                config_json=config_json, credentials_json=credentials_json
            )
        if cloud == "azure":
            self.dep = Azure(config_json=config_json, credentials_json=credentials_json)

    def deploy(self, model_id: str):
        """
        Deploy the model to the specified cloud platform.

        Parameters
        ----------
        model_id : str
            The ID of the model to be deployed.
        """
        self.dep.deploy(model_id)
