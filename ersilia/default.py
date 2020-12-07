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

# EOS conda
_resolve_script = "conda_env_resolve.py"
resolve_script = os.path.join(EOS, _resolve_script)
if not os.path.exists(resolve_script):
    shutil.copyfile(os.path.join(ROOT, "utils", "aux", _resolve_script), resolve_script)

snippet = """
# >>> ersilia >>>
eosconda() {
    EOS_MODEL_ENV=$(python %s $2);
    conda $1 $EOS_MODEL_ENV
}
# <<< ersilia <<<
""" % resolve_script


def bashrc_eosconda_snippet(overwrite=True):

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
        f.write(snippet.rstrip().lstrip())


if __name__ == "__main__":
    bashrc_eosconda_snippet()
