# TODO This is a placeholder for the environment variables module. This module will be able to manage environment variables, including secrets and API keys, and put them 
# in the right place for the model to use them. This module will also be able to get the environment variables from the right place, including the environment, the terminal, etc.
import os
import shutil
import subprocess
from ...db.environments.managers import DotenvManager
from ... import ErsiliaBase
from ...default import DOTENV_FILE


class GetEnvironmentVariables(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.logger.debug("Environmental variable getter created")

    @staticmethod
    def _parse_dotenv(dotenv_file):
        with open(dotenv_file, "r") as f:
            lines = f.readlines()
        env_dict = {}
        for line in lines:
            key, value = line.split("=")
            env_dict[key] = value
        return env_dict

    def _get_from_python(self, env):
        return os.environ.get(env)

    def _get_from_terminal(self, env):
        return os.system(f"echo ${env}")
    
    def _get_from_dotenv_in_cwd(self, env):
        dotenv_file = os.path.join(os.getcwd(), DOTENV_FILE)
        if os.path.exists(dotenv_file):
            env_dict = self.parse_dotenv(dotenv_file)
            if env in env_dict:
                return env_dict[env]
            else:
                return None
        else:
            return None
    
    def _get_from_dotenv_in_eos(self, env):
        return DotenvManager().get(env)
    
    def get_one_variable(self, env):
        self.logger.debug("Getting environmental variable: {0}".format(env))
        value = self._get_from_python(env)
        if value is not None:
            return value
        value = self._get_from_terminal(env)
        if value is not None:
            return value
        value = self._get_from_dotenv_in_cwd(env)
        if value is not None:
            return value
        value = self._get_from_github_in_eos(env)
        if value is not None:
            return value
        return None


# TODO: It is not clear this is what we need or want. We need to discuss this further since, in the case of Docker containers, we will have a risk of exposing secrets.
class PutEnvironmentVariables(ErsiliaBase):
    def __init__(self, model_id, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json, credentials_json=None)
        self.model_id = model_id
        self.logger.debug("Environmental variable setter created")

    def _get_framework_folder(self):
        bundle_path = self._get_bundle_location(model_id=self.model_id)
        bentoml_style_path = os.path.join(bundle_path, self.model_id, "artifact", "model", "framework")
        if os.path.exists(bentoml_style_path):
            return bentoml_style_path
        fastapi_style_path = os.path.join(bundle_path, self.model_id, "model", "framework")
        if os.path.exists(fastapi_style_path):
            return fastapi_style_path
        return None        

    def _send_to_framework_in_docker_as_dotenv(self, dotenv_file):
        # TODO: Implement this
        subprocess.run("", shell=True).wait()

    def _send_to_local_framework_as_dotenv(self, dotenv_file):
        # TODO: Test this
        framework_dir = self._get_framework_folder()
        if framework_dir is None:
            self.logger.debug("Framework directory not found. Cannot copy the .env file.")
            return None
        self.logger.debug("Copying .env file to the framework directory: {0}".format(framework_dir))
        shutil.copy(dotenv_file, os.path.join(framework_dir, DOTENV_FILE))

    def put(self, env_dict: dict):
        pass