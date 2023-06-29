Build container:

```bash
docker build -t ersilia .
```

Run in detached mode:

```bash
docker run -d --name ersilia_docker ersilia --entrypoint /bin/bash
```

To enter docker container, just run:

```bash
docker run -it --name ersilia_docker --entrypoint /bin/bash ersilia
```

You can also create an alias in your `.bashrc` (Linux) or `.bash_profile` (MacOSX):
