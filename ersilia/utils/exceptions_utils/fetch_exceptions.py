from .exceptions import ErsiliaError


class FetchErsiliaError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model: {0}".format(self.model_id)
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class InvalidUrlError(ErsiliaError):
    def __init__(self, url):
        self.message = "Provided URL is invalid: {0}".format(url)
        self.hints = "Open a browser and check that the URL is valid. You should see an API interface."
        ErsiliaError.__init__(self, self.message, self.hints)


class S3DownloaderError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model from S3: {0}".format(
            self.model_id
        )
        self.hints = (
            "Make sure that this model is actually in Ersilia Model Hub's S3 bucket"
        )
        ErsiliaError.__init__(self, self.message, self.hints)


class GetFetchErsiliaError(ErsiliaError):
    def __init__(self, model_id):
        self.model_id = model_id
        self.message = "Error occured while fetching model: {0}".format(self.model_id)
        self.hints = ""
        ErsiliaError.__init__(self, self.message, self.hints)


class FolderNotFoundError(ErsiliaError):
    def __init__(self, folder_name):
        self.folder_name = folder_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Folder " + self.folder_name + " not found!"
        return text

    def _get_hints(self):
        text = "Check that you can clone repositories from github\n"
        return text


class CondaEnvironmentExistsError(ErsiliaError):
    def __init__(self, env_name):
        self.environment_name = env_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Environment " + self.environment_name + " does not exist!"
        return text

    def _get_hints(self):
        text = "Try to look at the log files and see if there is an error in the conda installation\n"
        text += (
            "Please report the error at:\n - https://github.com/ersilia-os/ersilia\n"
        )
        return text


class ModelPackageInstallError(ErsiliaError):
    def __init__(self, package_name):
        self.package_name = package_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = (
            'Error occured while installing package by running "'
            + self.package_name
            + '" command \n'
        )
        return text

    def _get_hints(self):
        text = "Try to manually activate the model environment and install package manually.\n"
        text += "If this does not work, please report the error at:\n - https://github.com/ersilia-os/ersilia\n\n"
        return text


class VirtualEnvironmentSetupError(ErsiliaError):
    def __init__(self, venv_name):
        self.virtual_env_name = venv_name
        self.message = self._get_message()
        self.hints = self._get_hints()
        ErsiliaError.__init__(self, self.message, self.hints)

    def _get_message(self):
        text = "Virtual Environment " + self.virtual_env_name + " does not exist!"
        return text

    def _get_hints(self):
        text = "Try to look at the log files and see if there is an error in the venv installation\n"
        text += (
            "Please report the error at:\n - https://github.com/ersilia-os/ersilia\n"
        )
        return text
