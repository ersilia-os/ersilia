Build container:

```bash
docker build -t ersilia .
```

Run in detached mode:

```bash
docker run -d --name ersilia_docker ersilia
```

To enter docker container, just run:

```bash
docker exec -it ersilia_docker /bin/bash
```

You can also create an alias in your `.bashrc` (Linux) or `.bash_profile` (MacOSX):
