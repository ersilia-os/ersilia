import os
import json
import collections
import tempfile
import requests
from ..utils.terminal import run_command
from .. import ErsiliaBase
from ..auth.auth import Auth


class ReadmeParser(ErsiliaBase):

    def __init__(self, config_json):
        ErsiliaBase.__init__(self, config_json=config_json)

    def _raw_readme_url(self, model_id):
        url = "https://raw.githubusercontent.com/ersilia-os/{0}/master/README.md".format(model_id)
        return url

    def _gh_view(self, model_id):
        tmp_folder = tempfile.mkdtemp()
        tmp_file = os.path.join(tmp_folder, "view.md")
        cmd = "gh repo view {0}/{1} > {2}".format("ersilia-os", model_id, tmp_file)
        run_command(cmd, quiet=True)
        with open(tmp_file, "r") as f:
            text = f.read()
        return text

    def _title(self, lines):
        """Title is determined by the first main header in markdown"""
        for l in lines:
            if l[:2] == "# ":
                s = l[2:].strip()
                return s

    def _description(self, lines):
        """Description is what comes after the title and before the next header"""
        text = "\n".join(lines)
        return text.split("# ")[1].split("\n")[1].split("#")[0].strip()

    def _model_github_url(self, model_id):
        return "https://github.com/ersilia-os/{0}".format(model_id)

    def parse(self, model_id):
        readme = os.path.join(self._dest_dir, model_id, "README.md")
        if os.path.exists(readme):
            with open(readme, "r") as f:
                text = f.read()
        else:
            if Auth().is_contributor():
                text = self._gh_view(model_id)
                if not text:
                    return None
                text = "--".join(text.split("--")[1:])
            else:
                r = requests.get(self._raw_readme_url(model_id))
                if r.status_code != 200:
                    return None
                text = r.text
        lines = text.split(os.linesep)
        title = self._title(lines)
        descr = self._description(lines)
        results = {
            "model_id": model_id,
            "title": self._title(lines),
            "description": self._description(lines),
            "github_url": self._model_github_url(model_id)
        }
        return results


class ModelCard(object):

    def __init__(self, config_json=None):
        self.rp = ReadmeParser(config_json=config_json)

    def get(self, model_id, as_json=False):
        r = self.rp.parse(model_id)
        if as_json:
            return json.dumps(r, indent=4)
        else:
            return r
