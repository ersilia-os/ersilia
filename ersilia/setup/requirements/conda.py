from ...utils.terminal import run_command_check_output
from ... import throw_ersilia_exception
from ...utils.exceptions_utils.setup_exceptions import CondaSetupError

class CondaRequirement(object):
    def __init__(self):
        self.name = "conda"

    def is_installed(self):
        cmd = "command -v {0}".format(self.name)
        output = run_command_check_output(cmd)
        if self.name in output:
            return True
        else:
            raise CondaSetupError

    def install(self):
        pass
        # TODO
