"""See available models in the Ersilia Model Hub"""

from github import Github
from ..utils.config import Config


class List(object):

    def __init__(self, config=None):
        self.conf = Config(config)

    def spreadsheet(self):
        """List models available in our spreadsheets"""
        pass

    def github(self):
        """List models available in the GitHub model hub repository"""
        repo = Github().get_repo(self.conf.HUB.REPO)
        print(repo)

    def hub(self):
        """List models as available in our model hub repository"""
        pass

    def bento(self):
        """List models available as BentoServices"""
        pass

    def local(self):
        """List models as available in the local computer"""
        pass
