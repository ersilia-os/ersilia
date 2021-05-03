# Containers

Ersilia relies on Docker containers. Ersilia stores containers at [Docker Hub](https://hub.docker.com/orgs/ersiliaos).

## Upload image to the Ersilia Hub

### Prepare machine

```text
docker login --username=$USER --password=$PWD
```

### Write dockerfile

```text
FROM ersilia/basic

WORKDIR /app
ADD . /app
EXPOSE 80
```

### Build image

```text
docker build -t $DOCKER_ACC/$DOCKER_REPO:$IMG_TAG .
```

### Push image

```text
docker push $DOCKER_ACC/$DOCKER_REPO:$IMG_TAG
```

