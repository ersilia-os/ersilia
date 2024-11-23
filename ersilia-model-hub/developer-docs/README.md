---
description: >-
  In-depth documentation of the Ersilia Model Hub to support developers
  contribute to our open source platform
---

# Developer Docs

## Installation

To start contributing to Ersilia, please fork our main branch and work on it locally. We recommend installing Ersilia in editable mode:

```bash
conda create -n ersilia python=3.12
conda activate ersilia
git clone 
https://github.com/username/ersilia
 #your fork
pip install -e . 
```

Additional requirements:

* DockerHub: all models incorporated in Ersilia will be Dockerized for easy deployment. While the dockerization step happens in the cloud (GitHub Actions) it is recommended to install docker for model testing purposes.
* Git-LFS: many model’s checkpoints are too large (>100MB) for GitHub storage. Make sure that [git-lfs](https://git-lfs.com/) is installed and active in your system to push large files to the model repository.&#x20;
