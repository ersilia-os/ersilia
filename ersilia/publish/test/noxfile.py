import nox, os, sys
from pathlib import Path

ROOT = Path.home() / "eos"
TEMP = ROOT / "temp"
TEMP.mkdir(parents=True, exist_ok=True)
PWD = Path(__file__).resolve().parent.parent.parent
nox.options.envdir = str(TEMP / ".nox")

def run(s, *args):
  cmd = " ".join(str(a) for a in args)
  s.log(f"→ {cmd}")
  try:
    s.run(*args, external=True)
  except nox.command.CommandFailed as e:
    s.error(f"fail: {cmd}\n{e}")

def shim_bash(s):
  d = TEMP / "shims"
  d.mkdir(parents=True, exist_ok=True)
  p = d / "bash"
  p.write_text('#!/bin/sh\nexec /bin/sh "$@"\n')
  os.chmod(p, 0o755)
  s.env["PATH"] = str(d) + os.pathsep + s.env.get("PATH", "")


@nox.session(venv_backend='conda', python='3.12', reuse_venv=True)
def test_model_image(session):
    model_id = session.posargs[0] if len(session.posargs) >= 1 else "eos3b5e"
    version = session.posargs[1] if len(session.posargs) >= 2 else "dev"
    session.env["ERSILIA_NOX"] = "1"
    b = Path(sys.executable).parent
    session.env["PATH"] = str(b) + os.pathsep + session.env.get("PATH", "")
    session.env["PYTHON"] = sys.executable
    session.env["CONDA_PREFIX"] = str(Path(sys.prefix))
    shim_bash(session)
    print("PATH:", session.env["PATH"])
    print("PYTHON:", session.env["PYTHON"])
    print("CONDA_PREFIX:", session.env["CONDA_PREFIX"])
    run(
    session,
    "git",
    "clone",
    "--depth",
    "1",
    f"https://github.com/ersilia-os/ersilia.git",
    f"{TEMP}/ersilia"
  )
    session.install("-e", f"{TEMP}/ersilia[test]")
    if not model_id or not version:
        session.error("MODEL_ID and MODEL_VERSION must be set")

    run(
       session,
        "ersilia",
        "-v",
        "test",
        model_id,
        "--deep",
        "--from_dockerhub",
        "--version",
        version,
    )
