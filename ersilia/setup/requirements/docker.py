from ...utils.terminal import run_command_check_output
from ...utils.docker import is_inside_docker, resolve_platform


class DockerRequirement(object):
    def __init__(self):
        self.name = "docker"

    def is_inside_docker(self):
        return is_inside_docker()

    def is_installed(self):
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

    def install(self):
        # TODO
        msg = """
        Install docker
        """
        print(msg)
