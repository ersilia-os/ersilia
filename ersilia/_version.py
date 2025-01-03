import json
import os
from urllib.request import urlopen

STATIC_VERSION_FILE = "_static_version.py"
PACKAGE_NAME = "ersilia"

root = os.path.dirname(os.path.abspath(__file__))


def is_semver(tag):
    tag = tag.split(".")
    if len(tag) != 3:
        return False
    else:
        if tag[0].isnumeric() and tag[1].isnumeric() and tag[2].isnumeric():
            return True
        else:
            return False


def get_latest_semver_tag():
    owner = "ersilia-os"
    repo = "ersilia"
    with urlopen(f"https://api.github.com/repos/{owner}/{repo}/tags") as response:
        tags = json.load(response)
        if tags:
            for tag in tags:
                tag = tag["name"]
                if tag.startswith("v"):
                    tag = tag[1:]
                if is_semver(tag):
                    return tag
    return None


def increment_patch_version(version):
    version = version.split(".")
    version[2] = str(int(version[2]) + 1)
    return ".".join(version)


def get_version_for_setup():
    # version = get_latest_semver_tag()
    version = increment_patch_version(get_version_from_static())
    with open(os.path.join(root, STATIC_VERSION_FILE), "w") as f:
        f.write('version = "{0}"\n'.format(version))
    return version


def get_version_from_static():
    with open(os.path.join(root, STATIC_VERSION_FILE), "r") as f:
        return f.read().split('"')[1]


if os.path.exists(os.path.join(root, STATIC_VERSION_FILE)):
    __version__ = get_version_from_static()
else:
    __version__ = get_latest_semver_tag()
