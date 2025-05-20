import os
import platform
import stat
import tempfile

from .terminal import run_command


def is_inside_docker():
    """
    Check if the current environment is inside a Docker container.

    Returns
    -------
    bool
        True if inside a Docker container, False otherwise.
    """
    if os.path.isfile("/.dockerenv"):
        return True
    else:
        return False


class SystemChecker(object):
    """
    A class to check various system properties.

    Methods
    -------
    is_arm64()
        Check if the system architecture is ARM64.
    is_github_action()
        Check if the code is running inside a GitHub Action.
    is_inside_docker()
        Check if the code is running inside a Docker container.
    """

    def __init__(self):
        self.uname = platform.uname()

    def is_arm64(self):
        """
        Check if the system architecture is ARM64.

        Returns
        -------
        bool
            True if the system architecture is ARM64, False otherwise.
        """
        if self.uname.machine in ["arm64", "aarch64"]:
            return True
        if "arm64" in self.uname.version.lower():
            return True
        if "aarch64" in self.uname.version.lower():
            return True
        return False

    def is_github_action(self):
        """
        Check if the code is running inside a GitHub Action.

        Returns
        -------
        bool
            True if running inside a GitHub Action, False otherwise.
        """
        if os.environ.get("GITHUB_ACTIONS") == "true":
            return True
        else:
            return False

    def is_inside_docker(self):
        """
        Check if the code is running inside a Docker container.

        Returns
        -------
        bool
            True if running inside a Docker container, False otherwise.
        """
        return is_inside_docker()

    def get_country(self):
        """
        Get the country code of the system based on the IP address.

        Returns
        -------
        tuple
            The country of the system and its 3-character code (Country, Code).
        """

        bash_script = """#!/usr/bin/env bash

        # Usage: ./geoip.sh [IP_ADDRESS]
        IP="${1:-$(curl -s https://ifconfig.me)}"
        RESP=$(curl -s "http://ip-api.com/json/${IP}")

        # helper: extract the first JSON field value for a given key
        #   args: $1 = JSON text, $2 = field name
        get_field(){
        # strip up to "<key>":"..."
        local tmp="${1#*\\\"$2\\\":\\\"}"
        # then strip everything after the closing quote
        printf '%s' "${tmp%%\\\"*}"
        }

        status=$(get_field "$RESP" status)
        if [[ "$status" != success ]]; then
        msg=$(get_field "$RESP" message)
        echo "Lookup failed: ${msg:-$status}" >&2
        exit 1
        fi

        country=$(get_field "$RESP" country)
        code=$(get_field "$RESP" countryCode)

        echo "Country: $country"
        echo "Code:    $code"
        """

        tmp_dir = tempfile.mkdtemp(prefix="ersilia-")

        with open(os.path.join(tmp_dir, "geoip.sh"), "w") as f:
            f.write(bash_script)

        os.chmod(
            os.path.join(tmp_dir, "geoip.sh"),
            stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR,
        )

        run_command(f"bash {tmp_dir}/geoip.sh > {tmp_dir}/geoip.txt")

        with open(os.path.join(tmp_dir, "geoip.txt"), "r") as f:
            text = f.read()
            if "Country" in text:
                country = text.split("Country: ")[1].rstrip().split("\n")[0]
                code = text.split("Code:    ")[1].rstrip().split("\n")[0]
            else:
                country = "None"
                code = "None"

        return (country, code)
