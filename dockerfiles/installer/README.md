Build container:

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

You can also create an alias in your `.bashrc` (Linux) or `.bash_profile` (MacOSX):
