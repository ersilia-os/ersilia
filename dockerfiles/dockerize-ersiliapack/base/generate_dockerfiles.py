import os

CONDA_VERSIONS = [
    "py38_23.11.0-2",
    "py39_24.7.1-0", 
    "py310_24.7.1-0",
    "py311_24.7.1-0",
    "py312_24.7.1-0"
]

PIP_VERSIONS = [
    "3.8-slim-bullseye",
    "3.9-slim-bullseye",
    "3.10-slim-bullseye",
    "3.11-slim-bullseye",
    "3.12-slim-bullseye"
]

DOCKER_ENTRYPOINT = """
#!/bin/bash
set -ex
if [ -z "${MODEL}" ];
then
    echo "Model name has not been specified"
    exit 1
fi
ersilia_model_serve --bundle_path /root/bundles/$MODEL --port 3000
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

def generate_conda_dockerfile():
    conda_versions = CONDA_VERSIONS
    dockerfile = read_conda_base_dockerfile()    
    for version in conda_versions:
        version = version.strip()
        with open(os.path.join("Dockerfile.conda" + version), "w") as f:
            for line in dockerfile:
                f.write(line.replace("version", version))

def generate_pip_dockerfile():
    pip_versions = PIP_VERSIONS
    dockerfile = read_python_base_dockerfile()
    
    for version in pip_versions:
        version = version.strip()
        with open(os.path.join("Dockerfile.pip" + version), "w") as f:
            for line in dockerfile:
                f.write(line.replace("version", version))

if __name__ == "__main__":
    generate_conda_dockerfile()
    generate_pip_dockerfile()
    write_entrypoint()