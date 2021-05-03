import importools


class Rdkit(object):
    def __init__(self):
        try:
            self.rdkit = importlib.import_module("rdkit")
        except:
            pass

    def install(self):
        pass
