"""See available models in the Ersilia Model Hub"""

import subprocess
import pandas as pd
import webbrowser


class ModelList(object):

    def __init__(self):
        pass

    @staticmethod
    def spreadsheet():
        """List models available in our spreadsheets"""
        webbrowser.open("https://docs.google.com/spreadsheets/d/1WE-rKey0WAFktZ_ODNFLvHm2lPe27Xew02tS3EEwi28/edit#gid=1723939193") # TODO: do not just go to the website

    @staticmethod
    def github():
        """List models available in the GitHub model hub repository"""
        webbrowser.open("https://github.com/ersilia-os/") # TODO: do not just go to the website

    @staticmethod
    def hub():
        """List models as available in our model hub repository"""
        webbrowser.open("http://ersilia-os.github.io/ersilia-hub.github.io") # TODO: do not just go to the website

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

