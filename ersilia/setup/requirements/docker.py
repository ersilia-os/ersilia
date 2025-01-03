from ...utils.docker import resolve_platform
from ...utils.system import is_inside_docker
from ...utils.terminal import run_command_check_output


class DockerRequirement(object):
    """
    A class to handle the checking and management of Docker.

    Methods
    -------
    is_inside_docker()
        Checks if the current environment is inside a Docker container.
    is_working()
        Checks if Docker is working by running a test container.
    is_logged_in()
        Checks if Docker is logged in.
    is_active()
        Checks if Docker is active.
    is_installed()
        Checks if Docker is installed.
    """

    def __init__(self):
        self.name = "docker"

    def is_inside_docker(self) -> bool:
        """
        Checks if the current environment is inside a Docker container.

        Returns
        -------
        bool
            True if inside a Docker container, False otherwise.
        """
        return is_inside_docker()

    def is_working(self) -> bool:
        """
        Checks if Docker is working by running a test container.

        Returns
        -------
        bool
            True if Docker is working, False otherwise.
        """
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

    def is_logged_in(self) -> bool:
        """
        Checks if Docker is logged in.

        Returns
        -------
        bool
            True if Docker is logged in, False otherwise.
        """
        cmd = "docker info"
        output = run_command_check_output(cmd)
        if "Username" in output:
            return True
        else:
            return False

    def is_active(self) -> bool:
        """
        Checks if Docker is active.

        Returns
        -------
        bool
            True if Docker is active, False otherwise.
        """
        cmd = "docker ps"
        output = run_command_check_output(cmd)
        if "CONTAINER ID" in output:
            return True
        else:
            return False

    def is_installed(self) -> bool:
        """
        Checks if Docker is installed.

        Returns
        -------
        bool
            True if Docker is installed, False otherwise.
        """
        cmd = "docker --version"
        output = run_command_check_output(cmd)
        if "Docker version" in output:
            return True
        else:
            return False
