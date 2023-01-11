from ... import throw_ersilia_exception
from ...utils.exceptions_utils.setup_exceptions import PingError
import requests


class PingRequirement(object):
    def __init__(self):
        pass

    @throw_ersilia_exception
    def is_connected(self):
        url = "http://www.google.com/"
        try:
            _ = requests.get(url)
            return True
        except:
            raise PingError
