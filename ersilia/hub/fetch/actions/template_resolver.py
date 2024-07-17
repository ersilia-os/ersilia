import os
import http.client
import urllib.parse

from . import BaseAction

from ....default import GITHUB_ORG, ALLOWED_API_NAMES


class TemplateResolver(BaseAction):
    def __init__(self, model_id, repo_path=None, config_json=None):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.repo_path = repo_path

    def _check_file_in_repo(self, file_path):
        file_path = os.path.join(self.repo_path, file_path)
        if os.path.exists(file_path):
            return True
        else:
            return False

    def _check_file_in_github(self, file_path):
        url = "https://raw.githubusercontent.com/{0}/{1}/main/{2}".format(
            GITHUB_ORG, self.model_id, file_path
        )
        parsed_url = urllib.parse.urlparse(url)
        conn = http.client.HTTPSConnection(parsed_url.netloc)
        try:
            conn.request("HEAD", parsed_url.path)
            response = conn.getresponse()
            return response.status == 200
        except Exception as e:
            return False
        finally:
            conn.close()

    def _check_file(self, file_path):
        if self.repo_path is not None:
            return self._check_file_in_repo(file_path)
        else:
            return self._check_file_in_github(file_path)

    def is_fastapi(self):
        if not self._check_file("Dockerfile") and not self._check_file("install.yml"):
            return False
        has_sh = False
        for allowed_api in ALLOWED_API_NAMES:
            if self._check_file("model/framework/{0}.sh".format(allowed_api)):
                has_sh = True
        if not has_sh:
            return False
        return True

    def is_bentoml(self):
        if not self._check_file("pack.py"):
            return False
        if not self._check_file("Dockerfile"):
            return False
        if not self._check_file("src/service.py"):
            return False
        return True
