from .identifiers import FileIdentifier
import subprocess


class Conda(object):

    def __init__(self):
        self.fi = FileIdentifier()

    def activate(self, environtment_yml_file):
        name = self.fi.encode(environtment_yml_file)
        # TODO


"""
Create a tool that works like this:

ersilia-conda activate eos0aaa
# 1. Get checksum (obtained from environment.yml)
# 2. Check if it exists
# 3. Create conda environment if it doesn't exist
# 4. Save to ~/eos/conda_packages.db So that I can remove when I want...
# This will activate the conda environment
"""