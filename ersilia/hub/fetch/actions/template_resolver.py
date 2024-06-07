import os
import http.client
import urllib.parse

from . import BaseAction

from ....default import GITHUB_ORG


class TemplateResolver(BaseAction):
    def __init__(self, model_id, repo_path=None, config_json=None):
        BaseAction.__init__(
            self, model_id=model_id, config_json=config_json, credentials_json=None
        )
        self.repo_path = repo_path

    def _check_file_in_repo(self, file_path):
        if os.path.exists(os.path.join(self.repo_path, file_path)):
            return True
        else:
            return False

    def _check_file_in_github(self, file_path):
        url = "https://raw.githubusercontent.com/{0}/{1}/main/{2}".format(
            GITHUB_ORG, self.model_id, file_path
        )
        parsed_url = urllib.parse.urlparse(url)
        conn = http.client.HTTPConnection(parsed_url.netloc)
        try:
            conn.request("HEAD", parsed_url.path)
            response = conn.getresponse()
            return response.status == 200
        except Exception as e:
            print(f"Error checking URL: {e}")
            return False
        finally:
            conn.close()

    def _check_file(self, file_path):
        if self.repo_path is not None:
            return self._check_file_in_repo(file_path)
        else:
            return self._check_file_in_github(file_path)

    def is_fastapi(self):
        return self._check_file(
            "README.md"
        )  # TODO: For now, it is just a placeholder. It always returns True.

    def is_bentoml(self):
        return self._check_file("pack.yml")
