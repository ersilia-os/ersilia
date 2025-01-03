import importlib

from ...utils.terminal import run_command


class RdkitRequirement(object):
    """
    A class to handle the installation of the RDKit library.

    Methods
    -------
    install()
        Installs the RDKit library.
    """

    def __init__(self):
        self.name = "rdkit"
        try:
            importlib.import_module(self.name)
        except:
            self.install()

    def install(self) -> None:
        """
        Installs the RDKit library.

        This method installs the RDKit package using pip.

        Returns
        -------
        None
        """
        run_command("python -m pip install rdkit==2023.9.1")


class ChemblWebResourceClientRequirement(object):
    """
    A class to handle the installation of the ChEMBL web resource client library.

    Methods
    -------
    install()
        Installs the ChEMBL web resource client library.
    """

    def __init__(self):
        self.name = "chembl_webresource_client"
        try:
            importlib.import_module(self.name)
        except:
            self.install()

    def install(self) -> None:
        """
        Installs the ChEMBL web resource client library.

        This method installs the ChEMBL web resource client package using pip.

        Returns
        -------
        None
        """
        run_command("python -m pip install chembl_webresource_client")
