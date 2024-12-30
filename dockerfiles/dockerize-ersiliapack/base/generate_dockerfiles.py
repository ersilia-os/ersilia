import os
import sys

CONDA_VERSIONS = [
    "py38_23.11.0-2",
    "py39_24.7.1-0",
    "py310_24.7.1-0",
    "py311_24.7.1-0",
    "py312_24.7.1-0",
]

PIP_VERSIONS = [
    "3.8-slim-bullseye",
    "3.9-slim-bullseye",
    "3.10-slim-bullseye",
    "3.11-slim-bullseye",
    "3.12-slim-bullseye",
]

# By default, we generate all versions
version_to_build = sys.argv[1] if len(sys.argv) > 1 else "all"

# We serve the model at port 80 bec Ersilia port maps all containers to port 80
DOCKER_ENTRYPOINT = """#!/bin/bash
set -ex
if [ -z "${MODEL}" ];
then
    echo "Model name has not been specified"
    exit 1
fi
ersilia_model_serve --bundle_path /root/bundles/$MODEL --port 80
echo "Serving model $MODEL..."
"""


def read_conda_base_dockerfile():
    with open("Dockerfile.conda", "r") as f:
        return f.readlines()


def read_python_base_dockerfile():
    with open("Dockerfile.pip", "r") as f:
        return f.readlines()


def write_entrypoint():
    with open("docker-entrypoint.sh", "w") as f:
        f.write(DOCKER_ENTRYPOINT)


def generate_conda_dockerfile(version):
    dockerfile = read_conda_base_dockerfile()
    version = version.strip()
    with open(os.path.join("Dockerfile.conda" + version), "w") as f:
        for line in dockerfile:
            f.write(line.replace("version", version))


def generate_pip_dockerfile(version):
    dockerfile = read_python_base_dockerfile()
    version = version.strip()
    with open(os.path.join("Dockerfile.pip" + version), "w") as f:
        for line in dockerfile:
            f.write(line.replace("version", version))


if __name__ == "__main__":
    if version_to_build == "all":
        for version in CONDA_VERSIONS:
            generate_conda_dockerfile(version)
        for version in PIP_VERSIONS:
            generate_pip_dockerfile(version)
    else:
        if version_to_build in CONDA_VERSIONS:
            generate_conda_dockerfile(version_to_build)
        elif version_to_build in PIP_VERSIONS:
            generate_pip_dockerfile(version_to_build)
        else:
            print(
                "Invalid version specified. Please specify a valid version or 'all' to build all versions"
            )
    write_entrypoint()
