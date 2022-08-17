from ersilia.default import DEFAULT_DOCKER_PLATFORM
from ...utils.terminal import run_command_check_output
from ...utils.docker import is_inside_docker

from ...default import DEFAULT_DOCKER_PLATFORM


class DockerRequirement(object):
    def __init__(self):
        self.name = "docker"

    def is_inside_docker(self):
        return is_inside_docker()

    def is_installed(self):
        cmd = "docker run --platform {0} hello-world".format(DEFAULT_DOCKER_PLATFORM)
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
