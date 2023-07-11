Build container:

```bash
docker build -t ersiliaos/shell .
```

Or pull it:

```bash
docker pull ersiliaos/shell
```

Run in detached mode:

```bash
docker run -d --name ersilia_shell ersiliaos/shell
```

To enter docker container, just run:

```bash
docker run -it --name ersilia_docker --entrypoint /bin/bash ersilia
```
