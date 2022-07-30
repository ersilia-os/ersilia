"""See available models in the Ersilia Model Hub"""

import subprocess
import requests
import os
from .card import ModelCard
from ... import ErsiliaBase
from ...utils.identifiers.model import ModelIdentifier
from ...auth.auth import Auth
from ...default import GITHUB_ORG
from ... import logger

try:
    import webbrowser
except ModuleNotFoundError as err:
    webbrowser = None

try:
    from github import Github
except ModuleNotFoundError as err:
    Github = None

try:
    from tabulate import tabulate
except ModuleNotFoundError as err:
    tabulate = None


class CatalogTable(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def as_table(self):
        if not tabulate:
            return None
        else:
            return tabulate(
                self.data,
                headers=self.columns,
                tablefmt="fancy_grid",
                colalign=("center", "center", "center"),
            )

    def __str__(self):
        return self.as_table()

    def __repr__(self):
        return self.__str__()


class ModelCatalog(ErsiliaBase):
    def __init__(self, tabular_view=True, config_json=None):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.mi = ModelIdentifier()
        self.tabular_view = tabular_view

    def _is_eos(self, s):
        if self.mi.is_valid(s):
            if not self.mi.is_test(s):
                return True
        return False

    def _get_title(self, card):
        if "title" in card:
            return card["title"]
        if "Title" in card:
            return card["Title"]
        return None

    def _get_slug(self, card):
        if "slug" in card:
            return card["slug"]
        if "Slug" in card:
            return card["Slug"]
        return None

    def _get_mode(self, card):
        if "mode" in card:
            return card["mode"]
        if "Mode" in card:
            return card["Mode"]

    def airtable(self):
        """List models available in AirTable Ersilia Model Hub base"""
        if webbrowser:  # TODO: explore models online
            if not self.tabular_view:
                webbrowser.open(
                    "https://airtable.com/shr9sYjL70nnHOUrP/tblZGe2a2XeBxrEHP"
                )
            else:
                webbrowser.open("https://airtable.com/shrUcrUnd7jB9ChZV")

    def github(self):
        """List models available in the GitHub model hub repository"""
        if Github is None:
            token = None
        else:
            token = Auth().oauth_token()
        logger.debug(
            "Looking for model repositories in {0} organization".format(GITHUB_ORG)
        )
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
        logger.info("Found {0} models".format(len(models)))
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
            slug = self._get_slug(card)
            title = self._get_title(card)
            mode = self._get_mode(card)
            R += [[model_id, slug, title, mode]]
        return CatalogTable(R, columns=["MODEL_ID", "SLUG", "TITLE", "MODE"])

    def local(self):
        """List models available locally"""
        mc = ModelCard()
        mi = ModelIdentifier()
        R = []
        logger.debug("Looking for models in {0}".format(self._bundles_dir))
        for model_id in os.listdir(self._bundles_dir):
            if not self._is_eos(model_id):
                continue
            card = mc.get(model_id)
            slug = self._get_slug(card)
            title = self._get_title(card)
            mode = self._get_mode(card)
            R += [[model_id, slug, title, mode]]
        logger.info("Found {0} models".format(len(R)))
        return CatalogTable(data=R, columns=["MODEL_ID", "SLUG", "TITLE", "MODE"])

    def bentoml(self):
        """List models available as BentoServices"""
        result = subprocess.run(
            ["bentoml", "list"], stdout=subprocess.PIPE, env=os.environ
        )
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
            for i, idx in enumerate(zip(cut_idxs, cut_idxs[1:] + [None])):
                r += [row[idx[0] : idx[1]].rstrip()]
            R += [[r[0].split(":")[0]] + r]
        columns = ["MODEL_ID"] + columns
        return CatalogTable(data=R, columns=columns)
