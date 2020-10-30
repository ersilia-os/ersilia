"""See available models in the Ersilia Model Hub"""

import subprocess
import pandas as pd
import webbrowser
import os
from tabulate import tabulate
from .card import ModelCard
from .. import ErsiliaBase


class ModelCatalog(ErsiliaBase):

    def __init__(self, config_json=None, as_dataframe=True):
        ErsiliaBase.__init__(self, config_json=config_json)
        self.as_dataframe = as_dataframe

    def _return(self, df):
        if self.as_dataframe:
            return df
        h = list(df.columns)
        R = []
        for r in df.values:
            R += [r]
        return tabulate(R, headers=h)

    @staticmethod
    def spreadsheet():
        """List models available in our spreadsheets"""
        webbrowser.open("https://docs.google.com/spreadsheets/d/1WE-rKey0WAFktZ_ODNFLvHm2lPe27Xew02tS3EEwi28/edit#gid=1723939193") # TODO: do not just go to the website

    @staticmethod
    def github():
        """List models available in the GitHub model hub repository"""
        webbrowser.open("https://github.com/ersilia-os/") # TODO: do not just go to the website

    def hub(self):
        """List models available in our model hub repository"""
        mc = ModelCard()
        R = []
        for k, v in mc.data.items():
            R += [[k, v["title"], v["date"]]]
        df = pd.DataFrame(R, columns=["MODEL_ID", "TITLE", "DATE"])
        return self._return(df)

    def local(self):
        """List models available locally"""
        mc = ModelCard()
        R = []
        for model_id in os.listdir(self._bundles_dir):
            card = mc.get(model_id)
            R += [[model_id, card.title, card.date]]
        df = pd.DataFrame(R, columns=["MODEL_ID", "TITLE", "DATE"])
        return self._return(df)

    def bentoml(self):
        """List models available as BentoServices"""
        result = subprocess.run(['bentoml', 'list'], stdout=subprocess.PIPE)
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
            R += [r]
        df = pd.DataFrame(R, columns=columns)
        df["MODEL_ID"] = [x.split(":")[0] for x in list(df["BENTO_SERVICE"])]
        df = df[["MODEL_ID"]+columns]
        return self._return(df)

