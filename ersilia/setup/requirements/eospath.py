import os

from ... import throw_ersilia_exception
from ...default import EOS
from ...utils.exceptions_utils.setup_exceptions import EosHomePathNotFoundError


class EosHomePathRequirement(object):
    """
    A class to check the existence of the EOS home path.

    Methods
    -------
    eos_home_path_exists()
        Checks if the EOS home path exists.
    """

    def __init__(self):
        pass

    @throw_ersilia_exception()
    def eos_home_path_exists(self) -> bool:
        """
        Checks if the EOS home path exists.

        Returns
        -------
        bool
            True if the EOS home path exists, False otherwise.
        """
        if os.path.exists(EOS):
            return True
        else:
            raise EosHomePathNotFoundError
