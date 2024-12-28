from .exceptions import ErsiliaError

# ruff: noqa: D101, D102


class GitLfsSetupError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Git LFS is not installed! Git LFS is needed to download large files from our GitHub repositories."
        return text

    def _get_hints(self):
        text = "An easy way to install Git LFS is the following Conda command:\n"
        text += "$ conda install -c conda-forge git-lfs"
        return text


class GithubCliSetupError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "GitHub CLI is not installed! GitHub CLI is a fantastic tool to interact with GitHub. Ersilia uses it in the backend."
        return text

    def _get_hints(self):
        text = "An easy way to install the GitHub CLI is the following Conda command:\n"
        text += "$ conda install -c conda-forge gh"
        return text


class CondaSetupError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Conda is not installed! Conda is required to create virtual environments for each model."
        return text

    def _get_hints(self):
        text = "Visit https://docs.conda.io/projects/conda/en/latest/user-guide/install/ \n"
        text += "for information about installing Conda"
        return text


class PingError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "No internet connection. Internet connection is required for downloading models from GitHub repositories."
        return text

    def _get_hints(self):
        text = "Make sure that your computer is connected to the internet and try again. \n"
        return text


class EosHomePathNotFoundError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "EOS Home path not found. Looks like Ersilia is not installed correctly."
        return text

    def _get_hints(self):
        text = "Re-install Ersilia and try again. \n"
        return text
