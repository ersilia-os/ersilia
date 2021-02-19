# Containers

Ersilia relies on Docker containers.
Ersilia stores containers at [Docker Hub](https://hub.docker.com/orgs/ersiliaos).


## Upload image to the Ersilia Hub

### Prepare machine

```
docker login --username=$USER --password=$PWD
```

### Write dockerfile

```
FROM ersilia/basic

WORKDIR /app
ADD . /app
EXPOSE 80

```

### Build image
```
docker build -t $DOCKER_ACC/$DOCKER_REPO:$IMG_TAG .
```

### Push image

```
docker push $DOCKER_ACC/$DOCKER_REPO:$IMG_TAG
```
