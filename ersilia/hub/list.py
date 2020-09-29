"""See available models in the Ersilia Model Hub"""

import subprocess
import pandas as pd


class ModelList(object):

    def __init__(self):
        pass

    def spreadsheet(self):
        """List models available in our spreadsheets"""
        pass

    def github(self):
        """List models available in the GitHub model hub repository"""
        pass

    def hub(self):
        """List models as available in our model hub repository"""
        pass

    @staticmethod
    def bentoml():
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
        return df

    def local(self):
        """List models as available in the local computer"""
