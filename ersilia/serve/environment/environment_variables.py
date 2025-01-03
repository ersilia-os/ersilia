# TODO This is a placeholder for the environment variables module. This module will be able to manage environment variables, including secrets and API keys, and put them
# in the right place for the model to use them. This module will also be able to get the environment variables from the right place, including the environment, the terminal, etc.
import os
import shutil
import subprocess

from ... import ErsiliaBase

# from dotenv import load_dotenv
from ...default import DOTENV_FILE


class GetEnvironmentVariable(ErsiliaBase):
    """
    Class to get environment variables from various sources.

    Parameters
    ----------
    config_json : dict, optional
        Configuration in JSON format.
    """

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.logger.debug("Environmental variable getter created")

    def _get_from_python(self, env):
        return os.environ.get(env)

    def _get_from_terminal(self, env):
        return os.system(f"echo ${env}")

    def _get_from_dotenv_in_cwd(self, env):
        # load_dotenv()
        return os.environ.get(env)

    def _get_from_dotenv_in_eos(self, env):
        return os.environ.get(env)

    def _get_from_github_secrets(self, env):
        return

    def _get_from_bashrc(self, env):
        pass

    def get(self, env: str) -> str:
        """
        Get the value of an environment variable from various sources.

        Parameters
        ----------
        env : str
            The name of the environment variable.

        Returns
        -------
        str
            The value of the environment variable.
        """
        value = self._get_from_python(env)
        if value is not None:
            return value
        value = self._get_from_terminal(env)
        if value is not None:
            return value
        value = self._get_from_dotenv(env)
        if value is not None:
            return value
        value = self._get_from_github_secrets(env)
        if value is not None:
            return value
        value = self._get_from_bashrc(env)
        if value is not None:
            return value


class PutEnvironmentVariable(ErsiliaBase):
    """
    Class to set environment variables for a given model.

    Parameters
    ----------
    model_id : str
        The ID of the model.
    config_json : dict, optional
        Configuration in JSON format.
    """

    def __init__(self, model_id: str, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.logger.debug("Environmental variable setter created")

    def _get_framework_folder(self):
        bundle_path = self._get_bundle_location(model_id=self.model_id)
        bentoml_style_path = os.path.join(
            bundle_path, self.model_id, "artifact", "model", "framework"
        )
        if os.path.exists(bentoml_style_path):
            return bentoml_style_path
        fastapi_style_path = os.path.join(
            bundle_path, self.model_id, "model", "framework"
        )
        if os.path.exists(fastapi_style_path):
            return fastapi_style_path
        return None

    def _send_to_framework_in_docker_as_dotenv(self, dotenv_file):
        subprocess.run("", shell=True).wait()

    def _send_to_local_framework_as_dotenv(self, dotenv_file):
        framework_dir = self._get_framework_folder()
        if framework_dir is None:
            self.logger.debug(
                "Framework directory not found. Cannot copy the .env file."
            )
            return None
        self.logger.debug(
            "Copying .env file to the framework directory: {0}".format(framework_dir)
        )
        shutil.copy(dotenv_file, os.path.join(framework_dir, DOTENV_FILE))

    def put(self, env_dict: dict):
        """
        Set environment variables for the model.

        Parameters
        ----------
        env_dict : dict
            Dictionary of environment variables to set, where keys are variable names and values are variable values.
        """
        pass
