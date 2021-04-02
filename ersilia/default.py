from pathlib import Path
import shutil
import os

# EOS environmental variables
EOS = os.path.join(str(Path.home()), "eos")
if not os.path.exists(EOS):
    os.makedirs(EOS)
ROOT = os.path.dirname(os.path.realpath(__file__))
CHECKSUM_NCHAR = 8
CONDA_ENV_YML_FILE = "environment.yml"
DOCKERFILE_FILE = "Dockerfile"
GITHUB_ORG = "ersilia-os"
GITHUB_ERSILIA_REPO = "ersilia"
CONFIG_JSON = "config.json"
CREDENTIALS_JSON = "credentials.json"
INSTALL_STATUS_FILE = ".install.status"
DOCKER_BENTO_PATH = "/bento"
DOCKERHUB_ORG = "ersiliaos"
DOCKERHUB_LATEST_TAG = "latest"
DEFAULT_MODEL_ID = "eos0zzz"

# EOS conda
_resolve_script = "conda_env_resolve.py"
resolve_script = os.path.join(EOS, _resolve_script)
if not os.path.exists(resolve_script):
    shutil.copyfile(os.path.join(ROOT, "utils", "supp", _resolve_script), resolve_script)

snippet = """
# >>> ersilia >>>
# !! Contents within this block are managed by 'ersilia' !!
eosconda() {
    EOS_MODEL_ENV=$(python %s $2);
    conda $1 $EOS_MODEL_ENV
}

ersilia() {
    if [[ $1 == "conda" ]]; then
        eosconda "${@: 2}"
    elif [[ $1 == "auth" ]]; then
        gh auth "${@: 2}"
    else
        command ersilia "$@"
    fi
}
# <<< ersilia <<<
""" % resolve_script


def bashrc_cli_snippet(overwrite=True):
    """Write a conda snippet in the user profile.

    This function writes on the user profile to create an executable to work
    with conda environments based on model identifiers.

    Motivation behind this function is to define an ersilia CLI.

    Args:
        - overwrite (bool): Overwrite the current bash profile file if the eosconda string is found.
    """

    def bashrc_path():
        home_path = Path.home()
        rc = os.path.join(home_path, ".bashrc")
        if os.path.exists(rc):
            return rc
        pr = os.path.join(home_path, ".bash_profile")
        if os.path.exists(pr):
            return pr

    fn = bashrc_path()
    with open(fn, "r") as f:
        text = f.read()
    if snippet in text:
        if overwrite:
            text = text.split(snippet)[0] + text.split(snippet)[1]
        else:
            return
    with open(fn, "w") as f:
        f.write(text)
    with open(fn, "a+") as f:
        f.write(snippet)


if __name__ == "__main__":
    bashrc_cli_snippet()
