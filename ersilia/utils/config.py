"""Ersilia config.

The Config provide access to all sort of useful parameters.
"""
import os
import json
from ..default import (
    EOS,
    GITHUB_ORG,
    GITHUB_ERSILIA_REPO,
    CONFIG_JSON,
    CREDENTIALS_JSON,
)


SECRETS_JSON = "secrets.json"
GDRIVE_CLIENT_SECRETS_JSON = "gdrive_client_secrets.json"
ERSILIA_SECRETS_GITHUB_REPO = "ersilia-secrets"


class Checker(object):
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
        self._package_path()
        return self.development_path


class _Field(object):
    """Config Field placeholder."""

    def __init__(self, field_kv):
        """Initialize updating __dict__ and evaluating values."""
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
    """

    def __init__(self, json_file=None):
        """Initialize a Config instance.

        A Config instance is loaded from a JSON file.
        """
        if json_file is None:
            try:
                json_file = os.environ["EOS_CONFIG"]
            except KeyError as err:
                json_file = os.path.join(EOS, CONFIG_JSON)
            except Exception as err:
                raise err
        eval_obj_dict = _eval_obj(json_file)
        self.__dict__.update(eval_obj_dict)

    def keys(self):
        return self.__dict__.keys()


class Secrets(object):
    def __init__(self, overwrite=True):
        self.overwrite = overwrite
        self.secrets_json = os.path.join(EOS, SECRETS_JSON)
        self.gdrive_client_secrets_json = os.path.join(EOS, GDRIVE_CLIENT_SECRETS_JSON)

    def _fetch_from_github(self, remote_path, local_path):
        """Fetch filename from ersilia-secrets repository"""
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
        self._fetch_from_github(SECRETS_JSON, self.secrets_json)

    def fetch_gdrive_secrets_from_github(self):
        self._fetch_from_github(
            GDRIVE_CLIENT_SECRETS_JSON, self.gdrive_client_secrets_json
        )

    def to_credentials(self, json_file):
        """Convert secrets to credentials file"""
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
        # Local paths
        from .paths import Paths

        pt = Paths()
        local = {}
        # .. development models path
        dev_mod_path = pt.models_development_path()
        if dev_mod_path is None:
            v = "None"
        else:
            v = "'{0}'".format(dev_mod_path)
        local["DEVEL_MODELS_PATH"] = v
        cred["LOCAL"] = local
        with open(json_file, "w") as f:
            json.dump(cred, f, indent=4, sort_keys=True)
        return True


class Credentials(object):
    def __init__(self, json_file=None):
        if json_file is None:
            try:
                json_file = os.environ["EOS_CREDENTIALS"]
            except KeyError as err:
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
        return self.__dict__.keys()


__all__ = ["Config", "Credentials"]
