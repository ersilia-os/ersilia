from .exceptions import ErsiliaError


class GitLfsSetupError(ErsiliaError):
    def __init__(
        self,
    ):
        self.message = self._get_message()
        self.hints = self._get_hints()
        super().__init__(self.message, self.hints)

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
        super().__init__(self.message, self.hints)

    def _get_message(self):
        text = "GitHub CLI is not installed! GitHub CLI is a fantastic tool to interact with GitHub. Ersilia uses it in the backend."
        return text

    def _get_hints(self):
        text = "An easy way to install the GitHub CLI is the following Conda command:\n"
        text += "$ conda install -c conda-forge gh"
        return text


class CondaSetupError(ErsiliaError):
    pass


class PingError(ErsiliaError):
    pass


class EosHomePathNotFoundError(ErsiliaError):
    pass