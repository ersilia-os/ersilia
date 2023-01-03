from ...utils.terminal import run_command_check_output


class CondaRequirement(object):
    def __init__(self):
        self.name = "conda"

    def is_installed(self):
        cmd = "command -v {0}".format(self.name)
        output = run_command_check_output(cmd)
        if self.name in output:
            return True
        else:
            return False

    def install(self):
        pass
        # TODO
