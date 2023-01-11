from ... import throw_ersilia_exception
from ...default import EOS
from ...utils.exceptions_utils.setup_exceptions import EosHomePathNotFoundError
import os


class EosHomePathRequirement(object):
    def __init__(self):
        pass

    @throw_ersilia_exception
    def eos_home_path_exists(self):
        if os.path.exists(EOS):
            return True
        else:
            raise EosHomePathNotFoundError
