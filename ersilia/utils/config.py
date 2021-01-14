"""Ersilia config.

The Config provide access to all sort of useful parameters.
"""
import os
import json
from ..default import EOS, GITHUB_ORG
from autologging import logged
import requests

CONFIG_JSON = "config.json"
SECRETS_JSON = "secrets.json"
CREDENTIALS_JSON = "credentials.json"

ERSILIA_SECRETS_GITHUB_REPO = "ersilia-secrets"


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


@logged
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
                self.__log.debug("EOS_CONFIG environment variable not set. " + "Using default config file.")
                json_file = os.path.join(EOS, CONFIG_JSON)
            except Exception as err:
                raise err
        eval_obj_dict = _eval_obj(json_file)
        self.__dict__.update(eval_obj_dict)

    def keys(self):
        return self.__dict__.keys()


@logged
class Secrets(object):

    def __init__(self, overwrite=True):
        self.overwrite = overwrite
        self.secrets_json = os.path.join(EOS, SECRETS_JSON)

    def fetch_from_github(self):
        """Fetch secrets from ersilia-secrets repository"""
        from ..auth.auth import Auth
        auth = Auth()
        is_contributor = auth.is_contributor()
        if is_contributor:
            token = auth.oauth_token()
            from .download import GitHubDownloader
            ghd = GitHubDownloader(overwrite=self.overwrite, token=token)
            ghd.download_single(GITHUB_ORG, ERSILIA_SECRETS_GITHUB_REPO, SECRETS_JSON, self.secrets_json)

    def to_credentials(self, json_file):
        # TODO - Do something more complex than a mere copy.
        with open(self.secrets_json, "r") as f:
            sj = json.load(f)
        cred = {}
        for k,v in sj.items():
            cred[k] = "'{0}'".format(v)
        with open(json_file, "w") as f:
            json.dump(cred, f, indent=4)


@logged
class Credentials(object):

    def __init__(self, json_file=None):
        if json_file is None:
            try:
                json_file = os.environ["EOS_CREDENTIALS"]
            except KeyError as err:
                self.__log.debug("EOS_CREDENTIALS environment variable not set. " + "Using default credentials file.")
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


__all__ = [
    "Config",
    "Credentials"
]
