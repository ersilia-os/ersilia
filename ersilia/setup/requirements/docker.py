from ...utils.terminal import run_command_check_output
from ...utils.docker import resolve_platform
from ...utils.system import is_inside_docker


class DockerRequirement(object):
    def __init__(self):
        self.name = "docker"

    def is_inside_docker(self):
        return is_inside_docker()

    def is_working(self):
        cmd = "docker run --platform {0} --name my_hello_world hello-world".format(
            resolve_platform()
        )
        output = run_command_check_output(cmd)
        cmd = "docker rm -f my_hello_world"
        run_command_check_output(cmd)
        if "Hello from Docker!" in output:
            return True
        else:
            return False

    def is_logged_in(self):
        cmd = "docker info"
        output = run_command_check_output(cmd)
        if "Username" in output:
            return True
        else:
            return False

    def is_active(self):
        cmd = "docker ps"
        output = run_command_check_output(cmd)
        if "CONTAINER ID" in output:
            return True
        else:
            return False

    def is_installed(self):
        cmd = "docker --version"
        output = run_command_check_output(cmd)
        if "Docker version" in output:
            return True
        else:
            return False

    def is_logged_in(self):
        cmd = "docker info"
        output = run_command_check_output(cmd)
        if "Username" in output:
            return True
        else:
            return False
