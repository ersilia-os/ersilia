import json
import os

from ..default import (
    CONFIG_JSON,
    CREDENTIALS_JSON,
    EOS,
    GITHUB_ERSILIA_REPO,
    GITHUB_ORG,
)

SECRETS_JSON = "secrets.json"
GDRIVE_CLIENT_SECRETS_JSON = "gdrive_client_secrets.json"
ERSILIA_SECRETS_GITHUB_REPO = "ersilia-secrets"


class Checker(object):
    """
    A class to check and manage configuration and credentials files.

    Methods
    -------
    config()
        Ensure the configuration file exists.
    get_development_path()
        Get the development path.
    """

    def __init__(self):
        self.development_path = None
        self._config()
        self._credentials()

    def _package_path(self):
        if self.development_path is None:
            from .paths import Paths

            pt = Paths()
            self.development_path = pt.ersilia_development_path()

    def _config(self):
        dst = os.path.join(EOS, CONFIG_JSON)
        if os.path.exists(dst):
            return
        self._package_path()
        if self.development_path is None:
            src_exists = False
        else:
            src = os.path.join(self.development_path, CONFIG_JSON)
            src_exists = os.path.exists(src)
        if src_exists:
            os.symlink(src, dst)
        else:
            from .download import GitHubDownloader

            gd = GitHubDownloader(overwrite=True)
            gd.download_single(
                GITHUB_ORG,
                GITHUB_ERSILIA_REPO,
                CONFIG_JSON,
                os.path.join(EOS, CONFIG_JSON),
            )

    def _credentials(self):
        dst = os.path.join(EOS, CREDENTIALS_JSON)
        if os.path.exists(dst):
            return
        self._package_path()
        if self.development_path is None:
            src_exists = False
        else:
            src = os.path.join(self.development_path, CREDENTIALS_JSON)
            src_exists = os.path.exists(src)
        if src_exists:
            os.symlink(src, dst)
        else:
            sc = Secrets()
            sc.fetch_from_github()
            if self.development_path is None:
                done = sc.to_credentials(dst)
            else:
                done = sc.to_credentials(src)
                if done:
                    os.symlink(src, dst)

    def config(self):
        """
        Ensure the configuration file exists.

        If the configuration file does not exist, it will be created by copying from the development path or downloading from GitHub.
        """
        if os.path.exists(os.path.join(EOS, CONFIG_JSON)):
            return
        os.makedirs(EOS, exist_ok=True)
        self._package_path()
        dev_path = self.development_path
        if dev_path is not None:
            import shutil

            src = os.path.join(dev_path, CONFIG_JSON)
            dst = os.path.join(EOS, CONFIG_JSON)
            shutil.copyfile(src, dst)
        else:
            from .download import GitHubDownloader

            gd = GitHubDownloader(overwrite=True)
            gd.download_single(
                GITHUB_ORG,
                GITHUB_ERSILIA_REPO,
                CONFIG_JSON,
                os.path.join(EOS, CONFIG_JSON),
            )

    def get_development_path(self):
        """
        Get the development path.

        Returns
        -------
        str
            The development path.
        """
        self._package_path()
        return self.development_path


class _Field(object):
    def __init__(self, field_kv):
        tmp = dict()
        for k, v in field_kv.items():
            if type(v) == dict:
                tmp[k] = _Field(v)
            else:
                tmp[k] = eval(v)
        self.__dict__.update(tmp)

    def items(self):
        return self.__dict__.items()

    def asdict(self):
        return self.__dict__

    def __getitem__(self, key):
        return self.__dict__[key]


def _eval_obj(json_file):
    with open(json_file) as fh:
        obj_dict = json.load(fh)

    eval_obj_dict = dict()
    for k, v in obj_dict.items():
        if type(v) == dict:
            eval_obj_dict[k] = _Field(v)
        else:
            eval_obj_dict[k] = eval(v)
    return eval_obj_dict


class Config(object):
    """Config class.

    An instance of this object holds config file section as attributes.

    Parameters
    ----------
    json_file : str, optional
        The path to the JSON configuration file. Default is None.
    """

    def __init__(self, json_file=None):
        """Initialize a Config instance.

        A Config instance is loaded from a JSON file.
        """
        if json_file is None:
            try:
                json_file = os.environ["EOS_CONFIG"]
            except KeyError:
                json_file = os.path.join(EOS, CONFIG_JSON)
            except Exception as err:
                raise err
        eval_obj_dict = _eval_obj(json_file)
        self.__dict__.update(eval_obj_dict)

    def keys(self):
        """
        Get the keys of the configuration.

        Returns
        -------
        dict_keys
            The keys of the configuration.
        """
        return self.__dict__.keys()


class Secrets(object):
    """
    A class to manage secrets and credentials.

    Parameters
    ----------
    overwrite : bool, optional
        Whether to overwrite existing files. Default is True.
    """

    def __init__(self, overwrite=True):
        self.overwrite = overwrite
        self.secrets_json = os.path.join(EOS, SECRETS_JSON)
        self.gdrive_client_secrets_json = os.path.join(EOS, GDRIVE_CLIENT_SECRETS_JSON)

    def _fetch_from_github(self, remote_path, local_path):
        from ..auth.auth import Auth

        auth = Auth()
        is_contributor = auth.is_contributor()
        if is_contributor:
            token = auth.oauth_token()
            from .download import GitHubDownloader

            ghd = GitHubDownloader(overwrite=self.overwrite, token=token)
            ghd.download_single(
                GITHUB_ORG, ERSILIA_SECRETS_GITHUB_REPO, remote_path, local_path
            )

    def fetch_from_github(self):
        """
        Fetch secrets from the GitHub repository.
        """
        self._fetch_from_github(SECRETS_JSON, self.secrets_json)

    def fetch_gdrive_secrets_from_github(self):
        """
        Fetch Google Drive client secrets from the GitHub repository.
        """
        self._fetch_from_github(
            GDRIVE_CLIENT_SECRETS_JSON, self.gdrive_client_secrets_json
        )

    def to_credentials(self, json_file):
        """
        Convert secrets to a credentials file.

        Parameters
        ----------
        json_file : str
            The path to the credentials file.

        Returns
        -------
        bool
            True if the credentials file was created successfully, False otherwise.
        """
        if not os.path.exists(self.secrets_json):
            return False
        with open(self.secrets_json, "r") as f:
            sj = json.load(f)
        cred = {}
        # Start with secrets
        secrets = {}
        for k, v in sj.items():
            secrets[k] = "'{0}'".format(v)
        cred["SECRETS"] = secrets
        with open(json_file, "w") as f:
            json.dump(cred, f, indent=4, sort_keys=True)
        return True


class Credentials(object):
    """
    A class to manage credentials.

    Parameters
    ----------
    json_file : str, optional
        The path to the JSON credentials file. Default is None.
    """

    def __init__(self, json_file=None):
        if json_file is None:
            try:
                json_file = os.environ["EOS_CREDENTIALS"]
            except KeyError:
                json_file = os.path.join(EOS, CREDENTIALS_JSON)
            except Exception as err:
                raise err
        if os.path.exists(json_file):
            eval_obj_dict = _eval_obj(json_file)
            self.__dict__.update(eval_obj_dict)
            self.exists = True
        else:
            self.exists = False

    def keys(self):
        """
        Get the keys of the credentials.

        Returns
        -------
        dict_keys
            The keys of the credentials.
        """
        return self.__dict__.keys()


__all__ = ["Config", "Credentials"]
