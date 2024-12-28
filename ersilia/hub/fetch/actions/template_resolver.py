import http.client
import os
import urllib.parse

from ....default import ALLOWED_API_NAMES, GITHUB_ORG
from . import BaseAction


class TemplateResolver(BaseAction):
    """
    Resolves the template type (BentoML or FastAPI) for a model by checking the presence of specific files.

    Parameters
    ----------
    model_id : str
        Identifier of the model.
    repo_path : str, optional
        Path to the local repository, by default None.
    config_json : dict, optional
        Configuration settings for the resolver, by default None.

    Methods
    -------
    is_fastapi() -> bool
        Checks if the model uses or built with FastAPI.
    is_bentoml() -> bool
        Checks if the model uses or built with BentoML.
    """

    def __init__(self, model_id: str, repo_path: str = None, config_json: dict = None):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.repo_path = repo_path

    def _check_file_in_repo(self, file_path: str) -> bool:
        file_path = os.path.join(self.repo_path, file_path)
        if os.path.exists(file_path):
            return True
        else:
            return False

    def _check_file_in_github(self, file_path: str) -> bool:
        url = "https://raw.githubusercontent.com/{0}/{1}/main/{2}".format(
            GITHUB_ORG, self.model_id, file_path
        )
        parsed_url = urllib.parse.urlparse(url)
        conn = http.client.HTTPSConnection(parsed_url.netloc)
        try:
            conn.request("HEAD", parsed_url.path)
            response = conn.getresponse()
            return response.status == 200
        except Exception:
            return False
        finally:
            conn.close()

    def _check_file(self, file_path: str) -> bool:
        if self.repo_path is not None:
            return self._check_file_in_repo(file_path)
        else:
            return self._check_file_in_github(file_path)

    def is_fastapi(self) -> bool:
        """
        Checks if the model uses FastAPI.

        Returns
        -------
        bool
            True if the model uses FastAPI, False otherwise.
        """
        if not self._check_file("Dockerfile") and not self._check_file("install.yml"):
            return False
        has_sh = False
        for allowed_api in ALLOWED_API_NAMES:
            if self._check_file("model/framework/{0}.sh".format(allowed_api)):
                has_sh = True
        if not has_sh:
            return False
        return True

    def is_bentoml(self) -> bool:
        """
        Checks if the model uses BentoML.

        Returns
        -------
        bool
            True if the model uses BentoML, False otherwise.
        """
        if not self._check_file("pack.py"):
            return False
        if not self._check_file("Dockerfile"):
            return False
        if not self._check_file("src/service.py"):
            return False
        return True
