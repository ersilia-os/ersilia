# Copyright 2019 Atalaya Tech, Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import os
from collections import OrderedDict
from configparser import ConfigParser

from ..exceptions import BentoMLConfigException

logger = logging.getLogger(__name__)


class BentoMLConfigParser(ConfigParser):
    """
    BentoML configuration parser.

    Parameters
    ----------
    default_config : str
        Serve as default value when conf key not presented in environment var or user local config file.
    """

    def __init__(self, default_config: str, *args, **kwargs):
        ConfigParser.__init__(self, *args, **kwargs)

        if default_config is not None:
            self.read_string(default_config)

    @staticmethod
    def _env_var_name(section, key):
        return "BENTOML__{}__{}".format(section.upper(), key.upper())

    def get(self, section: str, key: str = None, **kwargs) -> str:  # pylint:disable=arguments-differ
        """
        A simple hierarchical config access, priority order:
        1. environment var
        2. user config file
        3. bentoml default config file

        Parameters
        ----------
        section : str
            The section of the configuration.
        key : str, optional
            The key within the section of the configuration. If not provided, the section is used as the key.

        Returns
        -------
        str
            The value of the configuration.

        Raises
        ------
        BentoMLConfigException
            If the section/key is not found in the configuration.
        """
        if key is None:
            key = section
            section = "core"
        section = str(section).lower()
        key = str(key).lower()

        env_var = self._env_var_name(section, key)
        if env_var in os.environ:
            return os.environ[env_var]

        if ConfigParser.has_option(self, section, key):
            return ConfigParser.get(self, section, key, **kwargs)
        else:
            raise BentoMLConfigException(
                "section/key '{}/{}' not found in BentoML config".format(section, key)
            )

    def as_dict(self, display_source: bool = False) -> dict:
        """
        Convert the configuration to a dictionary.

        Parameters
        ----------
        display_source : bool, optional
            If True, include the source of the configuration value (default is False).

        Returns
        -------
        dict
            The configuration as a dictionary.
        """
        cfg = {}

        for section in self:
            cfg.setdefault(section, OrderedDict())
            for k, val in self.items(section=section, raw=False):
                if display_source:
                    cfg[section][k] = (val, "<bentoml.cg>")
                else:
                    cfg[section][k] = val

        for ev in os.environ:
            if ev.startswith("BENTOML__"):
                _, section, key = ev.split("__")
                val = os.environ[ev]
                if display_source:
                    val = (val, "env var")
                cfg.setdefault(section.lower(), OrderedDict()).update(
                    {key.lower(): val}
                )

        return cfg

    def __repr__(self) -> str:
        """
        String representation of the BentoMLConfigParser.

        Returns
        -------
        str
            The string representation of the configuration.
        """
        return "<BentoML config: {}>".format(str(self.as_dict()))
