import requests

from ... import throw_ersilia_exception
from ...utils.exceptions_utils.setup_exceptions import PingError


class PingRequirement(object):
    """
    A class to check internet connectivity.

    Methods
    -------
    is_connected()
        Checks if the system is connected to the internet.
    """

    def __init__(self):
        pass

    @throw_ersilia_exception()
    def is_connected(self) -> bool:
        """
        Checks if the system is connected to the internet.

        Returns
        -------
        bool
            True if the system is connected to the internet, False otherwise.
        """
        url = "http://www.google.com/"
        try:
            _ = requests.get(url)
            return True
        except:
            raise PingError
