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

import json
from dataclasses import asdict, is_dataclass
from dataclasses import fields as get_fields


class DataclassJsonEncoder(json.JSONEncoder):
    """
    Special JSON encoder for dataclass types.

    Methods
    -------
    default(o)
        Override the default method to handle dataclass objects.
    """

    def default(self, o):  # pylint: disable=method-hidden
        """
        Override the default method to handle dataclass objects.

        Parameters
        ----------
        o : object
            The object to encode.

        Returns
        -------
        dict
            The encoded object as a dictionary.
        """
        if is_dataclass(o):
            if hasattr(o, "to_json"):
                return o.to_json()
            else:
                return asdict(o)
        return super().default(o)


class json_serializer:
    """
    A decorator for serializing dataclass objects to JSON.

    Parameters
    ----------
    fields : list, optional
        The fields to include in the JSON output (default is None, which includes all fields).
    compat : bool, optional
        If True, only include fields that have non-default values (default is False).

    Methods
    -------
    __call__(klass)
        Apply the decorator to the dataclass.
    """

    def __init__(self, fields=None, compat=False):
        self.fields = fields
        self.compat = compat

    @staticmethod
    def _extract_nested(obj):
        if hasattr(obj, "to_json"):
            return obj.to_json()
        return obj

    def __call__(self, klass):
        """
        Apply the decorator to the dataclass.

        Parameters
        ----------
        klass : type
            The dataclass type to decorate.

        Returns
        -------
        type
            The decorated dataclass type.

        Raises
        ------
        TypeError
            If the klass is not a dataclass.
        """
        if not is_dataclass(klass):
            raise TypeError(
                f"{self.__class__.__name__} only accepts dataclasses, "
                f"got {klass.__name__}"
            )
        default_map = {
            f.name: f.default_factory() if callable(f.default_factory) else f.default
            for f in get_fields(klass)
        }
        if self.fields is None:
            self.fields = tuple(k for k in default_map.keys() if not k.startswith("_"))

        if self.compat:

            def to_json(data_obj):
                return {
                    k: self._extract_nested(getattr(data_obj, k))
                    for k in self.fields
                    if default_map[k] != getattr(data_obj, k)
                }

        else:

            def to_json(data_obj):
                return {
                    k: self._extract_nested(getattr(data_obj, k)) for k in self.fields
                }

        klass.to_json = to_json
        return klass
