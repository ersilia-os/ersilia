"""See available models in the Ersilia Model Hub"""

import subprocess
import requests
import os
from .card import ModelCard
from .. import ErsiliaBase
from ..utils.paths import Paths
from ..auth.auth import Auth
from ..default import GITHUB_ORG

try:
    import webbrowser
except ModuleNotFoundError as err:
    webbrowser = None

try:
    from github import Github
except:
    Github = None

try:
    from tabulate import tabulate
except:
    tabulate = None


class CatalogTable(object):

    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def as_table(self):
        if not tabulate:
            return None
        else:
            return tabulate(self.data, headers=self.columns)

    def __str__(self):
        return self.as_table()

    def __repr__(self):
        return self.__str__()


class ModelCatalog(ErsiliaBase):

    def __init__(self, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.eos_regex = Paths._eos_regex()

    def _is_eos(self, s):
        if self.eos_regex.match(s):
            return True
        else:
            return False

    @staticmethod
    def backlog():
        """List models available in our spreadsheets for future incorporation"""
        if webbrowser:
            webbrowser.open("https://docs.google.com/spreadsheets/d/1WE-rKey0WAFktZ_ODNFLvHm2lPe27Xew02tS3EEwi28/edit#gid=1723939193") # TODO: do not just go to the website

    def github(self):
        """List models available in the GitHub model hub repository"""
        if Github is None:
            token = None
        else:
            token = Auth().oauth_token()
        if token:
            g = Github(token)
            repo_list = [i for i in g.get_user().get_repos()]
            repos = []
            for r in repo_list:
                owner, name = r.full_name.split("/")
                if owner != GITHUB_ORG:
                    continue
                repos += [name]
        else:
            repos = []
            url = "https://api.github.com/users/{0}/repos".format(GITHUB_ORG)
            results = requests.get(url).json()
            for r in results:
                repos += [r["name"]]
        models = []
        for repo in repos:
            if self._is_eos(repo):
                models += [repo]
        return models

    def hub(self):
        """List models available in Ersilia model hub repository"""
        mc = ModelCard()
        models = self.github()
        R = []
        for model_id in models:
            card = mc.get(model_id)
            if card is None:
                continue
            R += [[model_id, card["title"]]]
        return CatalogTable(R, columns=["MODEL_ID", "TITLE"])

    def local(self):
        """List models available locally"""
        mc = ModelCard()
        R = []
        for model_id in os.listdir(self._bundles_dir):
            card = mc.get(model_id)
            R += [[model_id, card["title"]]]
        return CatalogTable(data=R, columns=["MODEL_ID", "TITLE"])

    def bentoml(self):
        """List models available as BentoServices"""
        result = subprocess.run(['bentoml', 'list'], stdout=subprocess.PIPE, env=os.environ)
        result = [r for r in result.stdout.decode("utf-8").split("\n") if r]
        if len(result) == 1:
            return
        columns = ["BENTO_SERVICE", "AGE", "APIS", "ARTIFACTS"]
        header = result[0]
        values = result[1:]
        cut_idxs = []
        for col in columns:
            cut_idxs += [header.find(col)]
        R = []
        for row in values:
            r = []
            for i, idx in enumerate(zip(cut_idxs, cut_idxs[1:]+[None])):
                r += [row[idx[0]:idx[1]].rstrip()]
            R += [[r[0].split(":")[0]] + r]
        columns = ["MODEL_ID"] + columns
        return CatalogTable(data=R, columns=columns)
