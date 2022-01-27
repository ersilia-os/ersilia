from ...utils.terminal import run_command_check_output, run_command


class DockerRequirement(object):
    def __init__(self):
        self.name = "docker"

    def is_installed(self):
        cmd = "docker run hello-world"
        output = run_command_check_output(cmd)
        if "Hello from Docker!" in output:
            return True
        else:
            return False

    def install(self):
        # TODO
        msg = """
        Install docker
        """
        print(msg)
