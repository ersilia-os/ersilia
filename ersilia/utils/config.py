"""Ersilia config.

The Config provide access to all sort of useful parameters.
"""
import os
import json
from autologging import logged


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
                json_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../config.json')

            except Exception as err:
                raise err

        self.__log.debug('Loading config from: %s' % json_file)
        with open(json_file) as fh:
            obj_dict = json.load(fh)

        eval_obj_dict = dict()
        for k, v in obj_dict.items():
            if type(v) == dict:
                eval_obj_dict[k] = _Field(v)
            else:
                eval_obj_dict[k] = eval(v)
        self.__dict__.update(eval_obj_dict)

    def keys(self):
        return self.__dict__.keys()


__all__ = [
    "Config"
]
