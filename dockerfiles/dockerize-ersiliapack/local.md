## Introduction

Ersilia Pack dockerization is divided into two steps - building the base image that has ersilia pack, and then building the model image with this base image. There are two flavors of base images, viz., conda, and simple python. Model images follow the same pattern.

### Local build - base

1. To build the base images, make sure that the template Dockerfiles, ie, `Dockerfile.pip`, and `Dockerfile.conda` are in the same directory where you run `generate_dockerfile.py`.

2. This generates Dockerfiles for the python versions and conda versions specified in the script. 

3. This generates the following conda based and pip based Dockerfiles:
    - Dockerfile.condapy38_23.11.0-2
    - Dockerfile.condapy39_24.7.1-0
    - Dockerfile.condapy310_24.7.1-0
    - Dockerfile.condapy311_24.7.1-0
    - Dockerfile.condapy312_24.7.1-0
    - Dockerfile.pip3.8-slim-bullseye
    - Dockerfile.pip3.9-slim-bullseye
    - Dockerfile.pip3.10-slim-bullseye
    - Dockerfile.pip3.11-slim-bullseye
    - Dockerfile.pip3.12-slim-bullseye

And an entrypoint script called `docker-entrypoint.sh`.

4. Copy the desired Dockerfile and `docker-entrypoint.sh` file into ersilia-pack directory.

5. This step assumes you are inside the ersilia-pack directory and you have completed Step 4. Build the base docker image, eg for Python 3.12 conda image, as follows:

```
docker build -f Dockerfile.pip3.12-slim-bullseye -t ersiliaos/ersiliapack-py312:latest .
```

### Local build - model

1. cd to the drectory where the model is located (./model_folder)

2. If the model doesn't have any conda dependencies, copy `Dockerfile.pip`, otherwise copy `Dockerfile.conda`. 

3. Specify the correct base image based on the Python version required by the model on the first line of the DockerFile and the Model ID on the secondline, remove the .pip or .conda ending.l ID

4. To build the image, run:

```
docker build -t ersiliaos/eosxxxx:latest .
```