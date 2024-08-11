"""See available models in the Ersilia Model Hub"""

import subprocess
import requests
import shutil
import os
import json
import csv
from .card import ModelCard
from ... import ErsiliaBase
from ...utils.identifiers.model import ModelIdentifier
from ...auth.auth import Auth
from ...default import GITHUB_ORG, BENTOML_PATH, MODEL_SOURCE_FILE
from ... import logger

try:
    import webbrowser
except ModuleNotFoundError as err:
    webbrowser = None

try:
    from github import Github
except ModuleNotFoundError as err:
    Github = None


class CatalogTable(object):
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns

    def as_list_of_dicts(self):
        R = []
        for r in self.data:
            d = {}
            for i, c in enumerate(self.columns):
                d[c] = r[i]
            R += [d]
        return R

    def as_json(self):
        R = self.as_list_of_dicts()
        return json.dumps(R, indent=4)

    def write(self, file_name):
        with open(file_name, "w") as f:
            if file_name.endswith(".csv"):
                delimiter = ","
            elif file_name.endswith(".tsv"):
                delimiter = "\t"
            else:
                return None
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(self.columns)
            for r in self.data:
                writer.writerow(r)

    def __str__(self):
        return self.as_json()

    def __repr__(self):
        return self.__str__()


class ModelCatalog(ErsiliaBase):
    def __init__(self, config_json=None, only_identifier=True):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.mi = ModelIdentifier()
        self.only_identifier = only_identifier

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

    def _get_status(self, card):
        if "status" in card:
            return card["status"]
        if "Status" in card:
            return card["Status"]
        return None

    def _get_input(self, card):
        if "input" in card:
            return card["input"][0]
        if "Input" in card:
            return card["Input"][0]
        return None

    def _get_output(self, card):
        if "output" in card:
            return card["output"][0]
        if "Output" in card:
            return card["Output"][0]
        return None
    
    def _get_model_source(self, model_id):
        model_source_file = os.path.join(self._model_path(model_id), MODEL_SOURCE_FILE)
        if os.path.exists(model_source_file):
            with open(model_source_file) as f:
                return f.read().rstrip()
        else:
            return None
        
    def _get_service_class(self, card):
        if "service_class" in card:
            return card["service_class"]
        if "Service_class" in card:
            return card["Service_class"]
        return None

    def airtable(self):
        """List models available in AirTable Ersilia Model Hub base"""
        if webbrowser:
            webbrowser.open("https://airtable.com/shrUcrUnd7jB9ChZV")

    def _get_all_github_public_repos(self):
        url = "https://api.github.com/users/{0}/repos".format(GITHUB_ORG)
        while url:
            response = requests.get(url, params={"per_page": 100})
            response.raise_for_status()
            yield from response.json()
            if "next" in response.links:
                url = response.links["next"]["url"]  # get the next page
            else:
                break  # no more pages, stop the loop

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
            self.logger.debug("Token provided: ***")
            g = Github(token)
            repo_list = [i for i in g.get_user().get_repos()]
            repos = []
            for r in repo_list:
                owner, name = r.full_name.split("/")
                if owner != GITHUB_ORG:
                    continue
                repos += [name]
        else:
            self.logger.debug("Token not provided!")
            repos = []
            for repo in self._get_all_github_public_repos():
                repos += [repo["name"]]
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
        if self.only_identifier:
            R = []
            for model_id in models:
                R += [[model_id]]
            return CatalogTable(R, columns=["Identifier"])
        else:
            R = []
            for model_id in models:
                card = mc.get(model_id)
                if card is None:
                    continue
                slug = self._get_slug(card)
                title = self._get_title(card)
                R += [[model_id, slug, title]]
            return CatalogTable(R, columns=["Identifier", "Slug", "Title"])

    def local(self):
        """List models available locally"""
        mc = ModelCard()
        R = []
        logger.debug("Looking for models in {0}".format(self._bundles_dir))
        if self.only_identifier:
            R = []
            for model_id in os.listdir(self._bundles_dir):
                if not self._is_eos(model_id):
                    continue
                R += [[model_id]]
            columns = ["Identifier"]
        else:
            for model_id in os.listdir(self._bundles_dir):
                if not self._is_eos(model_id):
                    continue
                card = mc.get(model_id)
                slug = self._get_slug(card["card"])
                title = self._get_title(card["card"])
                status = self._get_status(card["card"])
                inputs = self._get_input(card["card"])
                output = self._get_output(card["card"])
                model_source = self._get_model_source(model_id)
                service_class = self._get_service_class(card)
                R += [[model_id, slug, title, status, inputs, output, model_source, service_class]]
            columns = [
                "Identifier",
                "Slug",
                "Title",
                "Status",
                "Input",
                "Output",
                "Model Source",
                "Service Class",
            ]
        logger.info("Found {0} models".format(len(R)))
        if len(R) == 0:
            return CatalogTable(data=[], columns=columns)
        return CatalogTable(data=R, columns=columns)

    def bentoml(self):
        """List models available as BentoServices"""
        try:
            result = subprocess.run(
                ["bentoml", "list"], stdout=subprocess.PIPE, env=os.environ, timeout=10
            )
        except Exception as e:
            shutil.rmtree(BENTOML_PATH)
            return None
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
        columns = ["Identifier"] + columns
        return CatalogTable(data=R, columns=columns)
