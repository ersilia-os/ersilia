---
description: Setting up Ersilia versions
---

# Versioning

## How to create a new version of Ersilia

### Version and tag semantics

Ersilia uses [semantic versioning](https://semver.org/), i.e. `{major}.{minor}.{patch}`. **GitHub tags and Ersilia versions are synchronized**. This means that Ersilia `0.1.14` corresponds tot the GitHub tag `v0.1.14`.

### GitHub tags

Ersilia versioning is **fully managed through GitHub**. Therefore, to create a new version of Ersilia, all you need to do is create a GitHub tag.

We recommend using **GitHub Desktop** to [tag a specific commit](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/managing-commits/managing-tags-in-github-desktop). Tags mush follow the `v0.0.0` notation.

Tags can also be created [from the commandline](https://git-scm.com/book/en/v2/Git-Basics-Tagging) as follows:

```bash
git tag v0.1.14
git push origin v0.1.14
```

### GitHub Actions workflows triggered

When a tag is pushed to the `master`  branch, the following actions will be triggered:

* [Upload Python package to PyPi](https://github.com/ersilia-os/ersilia/actions/workflows/python-publish.yml): This workflow will create a PyPi [Ersilia package](https://pypi.org/project/ersilia/) with version `0.1.14`.
* [Build and push Ersilia base image to DockerHub](https://github.com/ersilia-os/ersilia/blob/master/.github/workflows/ersilia-base-image-to-dockerhub.yml): This workflow will build and push an Ersilia [base image](https://hub.docker.com/repository/docker/ersiliaos/base/general) to DockerHub.



